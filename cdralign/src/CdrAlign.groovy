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
println "[CDRALIGN] Loading database"

def recordMap = new HashMap<String, Record>()

def antigenCountMap = new HashMap<String, Integer>()

def minDbScore = 2
def firstLine = true
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (it[-1].toInteger() >= minDbScore) {
            antigenCountMap.put(it[3], (antigenCountMap[it[3]] ?: 0) + 1)
        }
    } else {
        firstLine = false
    }
}

firstLine = true
def minCdr3CountPerAntigen = 10
new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (it[-1].toInteger() >= minDbScore && antigenCountMap[it[3]] >= minCdr3CountPerAntigen) {
            def record
            recordMap.put(it[1], record = (recordMap[it[0]] ?: new Record(it[0], it[1])))
            record.antigen.add(new AminoAcidSequence(it[3]))
        }
    } else {
        firstLine = false
    }
}

println "[CDRALIGN] Loaded ${recordMap.size()} unique CDR3s."

// align all-vs-all
def searchParameters = new TreeSearchParameters(7, 3, 3, 10)

def treeMap = new SequenceTreeMap<AminoAcidSequence, Record>(AminoAcidSequence.ALPHABET)

recordMap.values().each {
    treeMap.put(it.cdr3, it)
}

println "[CDRALIGN] Performing alignments."

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
            println "[CDRALIGN] Done all alignments for $counterVal records, " +
                    "total number of aligned CDR3 pairs is ${alignments.size()}"
        }
    }
}

println "[CDRALIGN] Done, ${alignments.size()} alignments performed. Proceeding to optimization"

// Run optimization
def listener = new ProgressListener() {
    @Override
    void progressUpdate(ProgressEvent event) {
        println "[MOEA stats]" + "\n" +
                "Time elapsed = " + event.elapsedTime + "\n" +
                "Number of function evaluations = " + event.currentNFE + "\n" +
                "Percent complete = " + event.percentComplete
    }
}
def result = new Executor()
        .distributeOnAllCores()
        .withProblem(new ScoringProblem(alignments))
        .withAlgorithm("NSGAII")
        .withProperty("populationSize", 200)
        .withMaxEvaluations(1000)
        .withProgressListener(listener).run()

//display the results
def alphabet = AminoAcidSequence.ALPHABET
def AAS = ['F', 'S', 'Y', 'C', 'W', 'L', 'P', 'H', 'Q', 'I',
           'M', 'T', 'N', 'K', 'R', 'V', 'A', 'D', 'E', 'G'] as List<Character>

new File("../solutions.txt").withPrintWriter { pw ->
    pw.println("id\tparameter\tfrom\tto\tvalue")
    result.eachWithIndex { Solution solution, int index ->
        def info = ScoringProblem.decode(solution)

        AAS.each { from ->
            AAS.each { to ->
                int score = info.scoring.getScore(alphabet.symbolToCode(from), alphabet.symbolToCode(to))
                pw.println(index + "\tsubstitution\t" + from + "\t" + to + "\t" + score)
            }
        }

        pw.println(index + "\tgap\tNA\tNA\t" + info.scoring.gapPenalty)
        pw.println(index + "\tmu\tNA\tNA\t" + info.mu)
        pw.println(index + "\tsigma\tNA\tNA\t" + info.sigma)
        pw.println(index + "\tthreshold\tNA\tNA\t" + info.threshold)
    }
}

new File("../roc.txt").withPrintWriter { pw ->
    pw.println("id\tsensitivity\tspecificity")
    result.eachWithIndex { Solution solution, int index ->
        pw.println(index + "\t" + solution.getObjective(0) + "\t" + solution.getObjective(1))
    }
}
