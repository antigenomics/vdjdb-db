import com.milaboratory.core.sequence.AminoAcidAlphabet
import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters

/**
 * Created by mikesh on 4/28/16.
 */

@Grapes(
        @Grab(group = 'com.milaboratory', module = 'milib', version = '1.3')
)

def species = "HomoSapiens", gene = "TRB"

def tcrAgPairSet = new HashSet<String>()

def notFirstLine = false

def agSet = new HashSet<String>(), agByTcrMap = new HashMap<String, String>()

new File("../database/vdjdb.txt").splitEachLine("\t") {
    if (notFirstLine) {
        if (it[1].equalsIgnoreCase(gene) && it[6].equalsIgnoreCase(species)) {
            tcrAgPairSet.add(it[2] + "\t" + it[10])
            agSet.add(it[10])
            agByTcrMap.put(it[2], it[10])
        }
    } else {
        notFirstLine = true
    }
}

new File("../net/").mkdir()

new File("../net/vdjdb.net.nodes.txt").withPrintWriter { pw ->
    pw.println("node\tag_type\tag_spec")
    agSet.sort().eachWithIndex { it, ind ->
        pw.println(it + "\t" + ind + "\t0")
    }
    agByTcrMap.each {
        pw.println(it.key + "\t0\t" + it.value)
    }
}

new File("../net/vdjdb.net.txt").withPrintWriter { pw ->
    pw.println("first_node\tsecond_node\tinteraction\tstrength")

    def tcrTree = new SequenceTreeMap<AminoAcidSequence, AminoAcidSequence>(AminoAcidSequence.ALPHABET)
    tcrAgPairSet.each {
        def (tcr, ag) = it.split("\t")

        pw.println(tcr + "\t" + ag + "\tta\t1.0")

        def tcrAA = new AminoAcidSequence(tcr)
        tcrTree.put(tcrAA, tcrAA)
    }

    def treeSearchParams = new TreeSearchParameters(2, 1, 1, 2)

    def processedEdges = new HashSet<String>()
    tcrTree.values().each { from ->
        def iter = tcrTree.getNeighborhoodIterator(from, treeSearchParams)

        def to
        while ((to = iter.next()) != null) {
            if (to != from) {
                def edge = from.toString() + "\t" + to.toString(),
                    edge_r = to.toString() + "\t" + from.toString(),
                    aln = iter.currentAlignment

                if (!(processedEdges.contains(edge) || processedEdges.contains(edge_r))) {
                    processedEdges.add(edge)
                    processedEdges.add(edge_r)

                    pw.println(edge + "\ttt\t" + aln.similarity())
                }
            }
        }
    }
}