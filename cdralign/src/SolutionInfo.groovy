import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

import static com.milaboratory.core.mutations.Mutation.getFrom;
import static com.milaboratory.core.mutations.Mutation.getTo;
import static com.milaboratory.core.mutations.Mutation.isDeletion;
import static com.milaboratory.core.mutations.Mutation.isInsertion;

public class SolutionInfo {
    private final LinearGapAlignmentScoring scoring;
    private final double threshold;

    SolutionInfo(LinearGapAlignmentScoring scoring, double threshold) {
        this.scoring = scoring;
        this.threshold = threshold;
    }

    public double computeScore(Alignment alignment) {
        Sequence reference = alignment.getSequence1();
        Mutations mutations = alignment.getAbsoluteMutations();
        double score = 0;

        for (int i = 0; i < reference.size(); i++) {
            byte aa = reference.codeAt(i);
            score += scoring.getScore(aa, aa);
        }

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i);

            if (isInsertion(mutation)) {
                score += scoring.getGapPenalty();
            } else {
                byte from = getFrom(mutation);
                score += isDeletion(mutation) ? scoring.getGapPenalty() :
                        (scoring.getScore(from, getTo(mutation)));
                score -= scoring.getScore(from, from);
            }
        }

        return score;
    }

    public LinearGapAlignmentScoring getScoring() {
        return scoring;
    }

    public double getThreshold() {
        return threshold;
    }
}
