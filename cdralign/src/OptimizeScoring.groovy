@Grapes(
        @Grab(group = 'org.moeaframework', module = 'moeaframework', version = '2.10')
)

import org.moeaframework.Executor
import org.moeaframework.core.Solution

import java.util.concurrent.ConcurrentLinkedQueue

def gene = args[0]

println "[OPTIMIZESCORING] Loading alignments"

def alignments = new ObjectInputStream(new FileInputStream("../" + gene + ".bin")).readObject() as ConcurrentLinkedQueue<RecordAlignment>

println "[OPTIMIZESCORING] Done. ${alignments.size()} alignments loaded."

def result = new Executor()
        .distributeOnAllCores()
        .withProblem(new ScoringProblem(alignments))
        .withAlgorithm("NSGAII")
        .withMaxEvaluations(10000)
        .run()

//display the results
System.out.format("Objective1  Objective2%n");

for (Solution solution : result) {
    System.out.format("%.4f      %.4f%n",
            solution.getObjective(0),
            solution.getObjective(1));
}