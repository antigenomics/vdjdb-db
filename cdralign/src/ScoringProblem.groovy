import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.alignment.LinearGapAlignmentScoring
import com.milaboratory.core.sequence.AminoAcidSequence
import org.moeaframework.core.Solution
import org.moeaframework.core.variable.EncodingUtils
import org.moeaframework.core.variable.RealVariable
import org.moeaframework.problem.AbstractProblem

import static com.milaboratory.core.mutations.Mutation.getFrom
import static com.milaboratory.core.mutations.Mutation.getTo
import static com.milaboratory.core.mutations.Mutation.isDeletion
import static com.milaboratory.core.mutations.Mutation.isInsertion

class ScoringProblem extends AbstractProblem {
    final Collection<RecordAlignment> alignments
    static final double MAX_DIAG = 1.0, MIN_NON_DIAG = -1.0, MIN_GAP = -1.0, VAR_FACTOR = 1000
    static final int N_SUBST = AminoAcidSequence.ALPHABET.size() * (AminoAcidSequence.ALPHABET.size() + 1) / 2,
                     N_SUBST_1 = AminoAcidSequence.ALPHABET.size(),
                     N_SUBST_2 = N_SUBST_1 * N_SUBST_1,
                     N_VARS = N_SUBST + 1

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
        def scoring = getScoring(solution)

        double overlapScore = 0.0, noOverlapScore = 0.0
        int nOverlap = 0

        alignments.each { RecordAlignment recordAlignment ->
            double score = computeScore(scoring, recordAlignment.alignment)

            if (recordAlignment.record1.antigen.any { recordAlignment.record2.antigen.contains(it) }) {
                overlapScore += score
                nOverlap++
            } else {
                noOverlapScore -= score
            }
        }

        solution.setObjective(0, overlapScore / nOverlap)
        solution.setObjective(1, noOverlapScore / (alignments.size() - nOverlap))
    }

    static LinearGapAlignmentScoring getScoring(Solution solution) {
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

        new LinearGapAlignmentScoring(AminoAcidSequence.ALPHABET, substitutionMatrix, gapPenalty)
    }

    static double computeScore(LinearGapAlignmentScoring scoring, Alignment alignment) {
        def mutations = alignment.absoluteMutations
        int score = 0
        for (int i = 0; i < mutations.size(); ++i) {
            int mutation = mutations.getMutation(i)

            if (isDeletion(mutation) || isInsertion(mutation))
                score += scoring.getGapPenalty()
            else //Substitution
                score += scoring.getScore(getFrom(mutation), getTo(mutation))
        }
        score / alignment.sequence1Range.length()
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

        solution
    }
}
