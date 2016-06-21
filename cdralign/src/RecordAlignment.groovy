import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

@Grapes(
        @Grab(group = 'com.milaboratory', module = 'milib', version = '1.3')
)

class RecordAlignment {
    final Record record1, record2
    final Alignment<AminoAcidSequence> alignment

    RecordAlignment(Record record1, Record record2, Alignment<AminoAcidSequence> alignment) {
        this.record1 = record1
        this.record2 = record2
        this.alignment = alignment
    }
}
