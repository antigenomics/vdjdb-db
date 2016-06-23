import com.milaboratory.core.sequence.AminoAcidSequence

@Grapes(
        @Grab(group = 'com.milaboratory', module = 'milib', version = '1.3')
)

class Record implements Serializable {
    final AminoAcidSequence cdr3
    final Gene gene
    final List<AminoAcidSequence> antigen = new ArrayList<>()

    Record(String gene, String cdr3) {
        this.gene = Gene."$gene"
        this.cdr3 = new AminoAcidSequence(cdr3)
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Record record = (Record) o

        if (cdr3 != record.cdr3) return false

        return true
    }

    int hashCode() {
        return cdr3.hashCode()
    }
}
