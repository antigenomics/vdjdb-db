@Grapes(
        [@Grab(group = 'com.milaboratory', module = 'milib', version = '1.3'),
                @Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.2.1'),
                @Grab(group = 'org.moeaframework', module = 'moeaframework', version = '2.10')]
)

import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool
import org.moeaframework.Executor
import org.moeaframework.core.Solution
import org.moeaframework.util.progress.ProgressEvent
import org.moeaframework.util.progress.ProgressListener

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

// requires a pre-built database
// load records
def sout = {
    println "[CDRALING ${new Date().toString()}] $it"
}

sout "Loading database"

def recordMap = new HashMap<String, Record>()

def antigenCountMap = new HashMap<String, Integer>()

def minDbScore = 2, species = "HomoSapiens"//, gene = "TRB"

def goodRecord = { List<String> record ->
    record[-1].toInteger() >= minDbScore && //record[0] == gene &&
            record[2] == species
}

def firstLine = true
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (goodRecord(it)) {
            antigenCountMap.put(it[3], (antigenCountMap[it[3]] ?: 0) + 1)
        }
    } else {
        firstLine = false
    }
}

firstLine = true
def minCdr3CountPerAntigen = 5
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (goodRecord(it) && antigenCountMap[it[3]] >= minCdr3CountPerAntigen) {
            def record
            recordMap.put(it[1], record = (recordMap[it[0]] ?: new Record(it[0], it[1])))
            record.antigen.add(new AminoAcidSequence(it[3]))
        }
    } else {
        firstLine = false
    }
}

sout "Loaded ${recordMap.size()} unique CDR3s."

// align all-vs-all
def searchParameters = new TreeSearchParameters(5, 2, 2, 7)

def treeMap = new SequenceTreeMap<AminoAcidSequence, Record>(AminoAcidSequence.ALPHABET)

recordMap.values().each {
    treeMap.put(it.cdr3, it)
}

sout "Performing alignments."

def alignments = new ConcurrentLinkedQueue<RecordAlignment>()
def counter = new AtomicInteger()

GParsPool.withPool(Runtime.getRuntime().availableProcessors()) {
    treeMap.values().eachParallel { Record from ->
        def iter = treeMap.getNeighborhoodIterator(from.cdr3, searchParameters)
        def to

        def toCdr3Hash = new HashSet<Record>()
        while ((to = iter.next()) != null) {
            if (from.cdr3 != to.cdr3 && from.gene == to.gene && !toCdr3Hash.contains(to)) {
                alignments.add(new RecordAlignment(from, to, iter.currentAlignment))
                toCdr3Hash.add(to)
            }
        }

        int counterVal
        if ((counterVal = counter.incrementAndGet()) % 100 == 0) {
            sout "Done all alignments for $counterVal records, " +
                    "total number of aligned CDR3 pairs is ${alignments.size()}"
        }
    }
}

sout "Done, ${alignments.size()} alignments performed. Proceeding to optimization"

// Run optimization
int popSize = 200, nGenerations = 1000
def listener = new ProgressListener() {
    @Override
    void progressUpdate(ProgressEvent event) {
        sout "[MOEA stats] Number of generations = " + (event.currentNFE / popSize)
    }
}

def problem = new ScoringProblem(alignments)

def result = new Executor()
        .distributeOnAllCores()
        .withProblem(problem)
        .withAlgorithm("NSGAII")
        .withProperty("populationSize", popSize)
        .withMaxEvaluations(nGenerations * popSize)
        .withProgressListener(listener).run()

//display the results
def alphabet = AminoAcidSequence.ALPHABET
def AAS = ['F', 'S', 'Y', 'C', 'W', 'L', 'P', 'H', 'Q', 'I',
           'M', 'T', 'N', 'K', 'R', 'V', 'A', 'D', 'E', 'G']

new File("../solutions.txt").withPrintWriter { pw ->
    pw.println("id\tparameter\tfrom\tto\tvalue")
    result.eachWithIndex { Solution solution, int index ->
        def info = problem.decode(solution)

        AAS.each { from ->
            AAS.each { to ->
                int score = info.scoring.getScore(alphabet.symbolToCode(from.charAt(0)),
                        alphabet.symbolToCode(to.charAt(0)))
                pw.println(index + "\tsubstitution\t" + from + "\t" + to + "\t" + score)
            }
        }

        pw.println(index + "\tgap\tNA\tNA\t" + info.scoring.gapPenalty)
        pw.println(index + "\tthreshold\tNA\tNA\t" + info.threshold)
        info.positionWeights.eachWithIndex { double value, int i ->
            pw.println(index + "\tposition_weight\t" + i + "\tNA\t" + value)
        }
    }
}

new File("../roc.txt").withPrintWriter { pw ->
    pw.println("id\tsensitivity\tspecificity")
    result.eachWithIndex { Solution solution, int index ->
        pw.println(index + "\t" + solution.getObjective(0) + "\t" + solution.getObjective(1))
    }
}
