import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.moeaframework.core.Solution
import org.moeaframework.core.variable.EncodingUtils
import org.moeaframework.core.variable.RealVariable
import org.moeaframework.problem.AbstractProblem

import static com.milaboratory.core.mutations.Mutation.getFrom
import static com.milaboratory.core.mutations.Mutation.getPosition
import static com.milaboratory.core.mutations.Mutation.getTo
import static com.milaboratory.core.mutations.Mutation.isDeletion
import static com.milaboratory.core.mutations.Mutation.isInsertion

class ScoringProblem extends AbstractProblem {
    final Collection<RecordAlignment> alignments
    static final double MAX_DIAG = 1.0, MIN_NON_DIAG = -1.0, MIN_GAP = -1.0, VAR_FACTOR = 1000,
                        SCORE_RANGE = 30 * Math.max(MAX_DIAG, Math.max(-MIN_NON_DIAG, -MIN_GAP))
    static final int N_SUBST = AminoAcidSequence.ALPHABET.size() * (AminoAcidSequence.ALPHABET.size() + 1) / 2,
                     N_SUBST_1 = AminoAcidSequence.ALPHABET.size(),
                     N_SUBST_2 = N_SUBST_1 * N_SUBST_1,
                     N_VARS = N_SUBST + 2

    ScoringProblem(Collection<RecordAlignment> alignments) {
        super(N_VARS, 2)
        this.alignments = alignments
    }

    /*
    @Override
    Solution generate() {
        Solution solution = newSolution()

        (0..<N_VARS).each { int i ->
            def var = solution.getVariable(i) as RealVariable
            var.setValue(PRNG.nextDouble(var.lowerBound, var.upperBound))
        }

        evaluate(solution)

        solution
    }*/

    @Override
    void evaluate(Solution solution) {
        def solutionInfo = decode(solution)

        int TP = 0, FP = 0, TN = 0, FN = 0

        alignments.each { RecordAlignment recordAlignment ->
            double score = solutionInfo.computeScore(recordAlignment.record1.cdr3,
                    recordAlignment.alignment)

            if (recordAlignment.record1.antigen.any { recordAlignment.record2.antigen.contains(it) }) {
                if (score >= solutionInfo.threshold) {
                    TP++
                } else {
                    FN++
                }
            } else {
                if (score >= solutionInfo.threshold) {
                    FP++
                } else {
                    TN++
                }
            }
        }

        solution.setObjective(0, -TP / (double) Math.max(1, TP + FN))
        solution.setObjective(1, -TN / (double) Math.max(1, TN + FP))
    }

    static SolutionInfo decode(Solution solution) {
        double[] vars = EncodingUtils.getReal(solution)

        int[] substitutionMatrix = new int[N_SUBST_2]

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                int var = VAR_FACTOR * vars[k]
                substitutionMatrix[i * N_SUBST_1 + j] = var
                if (i != j)
                    substitutionMatrix[j * N_SUBST_1 + i] = var
                k++
            }
        }

        int gapPenalty = VAR_FACTOR * vars[k]
        def scoring = new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET, substitutionMatrix, gapPenalty)

        def mu = vars[++k]
        def sigma = vars[++k]
        def threshold = VAR_FACTOR * vars[++k]

        new SolutionInfo(scoring, mu, sigma, threshold)
    }

    @Override
    Solution newSolution() {
        Solution solution = new Solution(N_VARS, 2)

        int k = 0
        for (int i = 0; i < N_SUBST_1; i++) {
            for (int j = i; j < N_SUBST_1; j++) {
                solution.setVariable(k, i == j ? new RealVariable(0, MAX_DIAG) :
                        new RealVariable(MIN_NON_DIAG, 0))
                k++
            }
        }

        solution.setVariable(k, new RealVariable(MIN_GAP, -2 / VAR_FACTOR))

        solution.setVariable(++k, new RealVariable(-3.0, 3.0))
        solution.setVariable(++k, new RealVariable(1.0, 10.0))
        solution.setVariable(++k, new RealVariable(-SCORE_RANGE, SCORE_RANGE))

        solution
    }
}
