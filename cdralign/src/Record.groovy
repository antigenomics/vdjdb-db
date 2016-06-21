import com.milaboratory.core.sequence.AminoAcidSequence

@Grapes(
        @Grab(group = 'com.milaboratory', module = 'milib', version = '1.3')
)

class Record {
    final AminoAcidSequence cdr3
    final List<AminoAcidSequence> antigen = new ArrayList<>()

    Record(String cdr3) {
        this.cdr3 = new AminoAcidSequence(cdr3)
    }
}
