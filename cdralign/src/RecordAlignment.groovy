import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class RecordAlignment implements Serializable {
    final Record record1, record2
    final Alignment<AminoAcidSequence> alignment
    final boolean antigensMatch

    RecordAlignment(Record record1, Record record2, Alignment<AminoAcidSequence> alignment) {
        this.record1 = record1
        this.record2 = record2
        this.alignment = alignment
        this.antigensMatch = record1.antigen.any { record2.antigen.contains(it) }
    }
}
