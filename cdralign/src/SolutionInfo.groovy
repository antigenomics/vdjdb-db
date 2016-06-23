import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.moeaframework.core.Solution
import org.moeaframework.core.variable.EncodingUtils

import static com.milaboratory.core.mutations.Mutation.getFrom
import static com.milaboratory.core.mutations.Mutation.getPosition
import static com.milaboratory.core.mutations.Mutation.getTo
import static com.milaboratory.core.mutations.Mutation.isDeletion
import static com.milaboratory.core.mutations.Mutation.isInsertion

class SolutionInfo {
    final LinearGapAlignmentScoring scoring
    final double mu, sigma, threshold

    SolutionInfo(LinearGapAlignmentScoring scoring, double mu, double sigma, double threshold) {
        this.scoring = scoring
        this.mu = mu
        this.sigma = sigma
        this.threshold = threshold
    }

    double computeScore(AminoAcidSequence reference,
                        Alignment alignment) {
        def mutations = alignment.absoluteMutations
        double score = 0

        double halfLength = reference.size() / 2.0

        for (int i = 0; i < reference.size(); i++) {
            byte aa = reference.codeAt(i)
            score += scoring.getScore(aa, aa) * computeWeight(i - halfLength)
        }

        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i)
            double weight = computeWeight(getPosition(mutation) - halfLength)

            double currentScore
            if (isInsertion(mutation)) {
                currentScore = scoring.getGapPenalty()
            } else {
                byte from = getFrom(mutation)
                currentScore = isDeletion(mutation) ? scoring.getGapPenalty() :
                        (scoring.getScore(from, getTo(mutation)))
                currentScore -= scoring.getScore(from, from)
            }

            score += currentScore * weight
        }

        score
    }

    static final double TWO_PI_SQRT = Math.sqrt(2 * Math.PI)

    double computeWeight(double x) {
        double delta = x - mu
        Math.exp(-delta * delta / sigma / sigma) / TWO_PI_SQRT / sigma
    }
}
