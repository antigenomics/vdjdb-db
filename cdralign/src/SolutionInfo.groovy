import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

import static com.milaboratory.core.mutations.Mutation.getFrom;
import static com.milaboratory.core.mutations.Mutation.getTo;
import static com.milaboratory.core.mutations.Mutation.isDeletion;
import static com.milaboratory.core.mutations.Mutation.isInsertion;

class SolutionInfo {
    final LinearGapAlignmentScoring scoring
    final double threshold

    SolutionInfo(LinearGapAlignmentScoring scoring, double threshold) {
        this.scoring = scoring
        this.threshold = threshold
    }

    public double computeScore(Alignment alignment) {
        Sequence reference = alignment.sequence1
        Mutations mutations = alignment.absoluteMutations
        double score = 0

        for (int i = 0; i < reference.size(); i++) {
            byte aa = reference.codeAt(i)
            score += scoring.getScore(aa, aa)
        }

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i)

            if (isInsertion(mutation)) {
                score += scoring.getGapPenalty()
            } else {
                byte from = getFrom(mutation)
                score += isDeletion(mutation) ? scoring.gapPenalty :
                        (scoring.getScore(from, getTo(mutation)))
                score -= scoring.getScore(from, from)
            }
        }

        score
    }
}
