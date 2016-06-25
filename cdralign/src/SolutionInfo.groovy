import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.alignment.LinearGapAlignmentScoring;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.Sequence;

import static com.milaboratory.core.mutations.Mutation.getFrom
import static com.milaboratory.core.mutations.Mutation.getPosition;
import static com.milaboratory.core.mutations.Mutation.getTo;
import static com.milaboratory.core.mutations.Mutation.isDeletion;
import static com.milaboratory.core.mutations.Mutation.isInsertion;

class SolutionInfo {
    final LinearGapAlignmentScoring scoring
    final double[] positionWeights
    final double threshold
    final int posCenterBin

    SolutionInfo(LinearGapAlignmentScoring scoring,
                 double[] positionWeights,
                 double threshold) {
        this.scoring = scoring
        this.positionWeights = positionWeights
        this.posCenterBin = positionWeights.length / 2
        this.threshold = threshold
    }


    private int getBin(int centeredPos) {
        int k = centeredPos + posCenterBin
        k < 0 ? 0 : (k < positionWeights.length ? k : (positionWeights.length - 1))
    }

    private double getPositionWeight(int pos, int cdr3Length) {
        int center = cdr3Length / 2
        if (cdr3Length % 2 == 0) {
            return 0.5 * positionWeights[getBin(pos - center)] + 0.5 * positionWeights[getBin(pos - center + 1)]
        } else {
            return positionWeights[getBin(pos - center)]
        }
    }

    double computeScore(Alignment alignment) {
        Sequence reference = alignment.sequence1
        Mutations mutations = alignment.absoluteMutations
        double score = 0

        for (int i = 0; i < reference.size(); i++) {
            byte aa = reference.codeAt(i)
            score += scoring.getScore(aa, aa) * getPositionWeight(i, reference.size())
        }

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i)

            double deltaScore = 0
            if (isInsertion(mutation)) {
                deltaScore += scoring.getGapPenalty()
            } else {
                byte from = getFrom(mutation)
                deltaScore += isDeletion(mutation) ? scoring.gapPenalty :
                        (scoring.getScore(from, getTo(mutation)))
                deltaScore -= scoring.getScore(from, from)
            }

            score += deltaScore * getPositionWeight(getPosition(mutation), reference.size())
        }

        score
    }
}
