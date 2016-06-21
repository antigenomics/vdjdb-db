import com.milaboratory.core.alignment.Alignment
@Grapes(
        [@Grab(group = 'com.milaboratory', module = 'milib', version = '1.3'),
                @Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.2.1')]
)



import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.atomic.AtomicInteger

// requires a pre-built database
// load records
println "[CDRALIGN] Loading database"

def gene = "TRB"
def recordMap = new HashMap<String, Record>()
def firstLine = true

new File("../../database/vdjdb.slim.txt").splitEachLine("\t") {
    if (!firstLine) {
        if (it[0] == gene) {
            def record
            recordMap.put(it[1], record = (recordMap[it[0]] ?: new Record(it[1])))
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

def chooseAlignment = { Alignment<AminoAcidSequence> a1, Alignment<AminoAcidSequence> a2 ->
    a1 ? a2 : (a1.score < a2.score ? a2 : a1)
}

GParsPool.withPool(Runtime.getRuntime().availableProcessors()) {
    treeMap.values().eachParallel { Record from ->
        def iter = treeMap.getNeighborhoodIterator(from.cdr3, searchParameters)
        def to
        def alignmentMap = new HashMap<Record, Alignment<AminoAcidSequence>>()

        while ((to = iter.next()) != null) {
            if (from.cdr3 != to.cdr3) {
                alignmentMap.put(to, chooseAlignment(alignmentMap[to], iter.currentAlignment))
            }
        }

        alignmentMap.each {
            alignments.add(new RecordAlignment(from, it.key, it.value))
        }

        int counterVal
        if ((counterVal = counter.incrementAndGet()) % 100 == 0) {
            println "[CDRALIGN] Done all alignments for $counterVal records, " +
                    "total number of aligned CDR3 pairs is ${alignments.size()}"
        }
    }
}

println "[CDRALIGN] Done, ${alignments.size()} alignments performed."

