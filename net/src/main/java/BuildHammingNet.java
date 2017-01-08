import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Stream;

import com.milaboratory.core.alignment.Alignment;
import com.milaboratory.core.mutations.MutationType;
import com.milaboratory.core.mutations.Mutations;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;
import com.milaboratory.core.tree.TreeSearchParameters;

public class BuildHammingNet {
    public static void main(String[] args) throws IOException {
        int maxSubstitutions = Integer.parseInt(args[0]),
                maxIndels = Integer.parseInt(args[1]);

        String inputFileName = args[2], outputFileName = args[3];

        final Map<AminoAcidSequence, Cdr3Info> cdr3AntigenMap = new HashMap<>();
        try (Stream<String> stream = Files.lines(new File(inputFileName).toPath())) {
            final boolean[] headerLine = {true};
            stream.forEach(line -> {
                if (headerLine[0]) {
                    headerLine[0] = false;
                } else {
                    String[] splitString = line.split("\t");
                    cdr3AntigenMap.compute(new AminoAcidSequence(splitString[0]),
                            (aminoAcidSequence, cdr3Info) -> {
                                if (cdr3Info == null) {
                                    cdr3Info = new Cdr3Info(aminoAcidSequence);
                                }
                                cdr3Info.addAntigen(splitString[1]);
                                return cdr3Info;
                            });
                }
            });
        }

        System.out.println("Loaded " + cdr3AntigenMap.size() + " cdr3 sequences.");

        final SequenceTreeMap<AminoAcidSequence, Cdr3Info> stm = new SequenceTreeMap<>(AminoAcidSequence.ALPHABET);

        cdr3AntigenMap.entrySet().forEach(kvp -> stm.put(kvp.getKey(), kvp.getValue()));

        final TreeSearchParameters tsp = new TreeSearchParameters(maxSubstitutions, maxIndels, maxIndels);

        String header = "cdr3.1\tcdr3.2\tsame.ag\tsubst\tins\tdel\tedge.id\tweight";

        try (final PrintWriter pw = new PrintWriter(new File(outputFileName))) {
            pw.println(header);

            final Queue<String> lines = new ConcurrentLinkedQueue<>();

            Thread writeThread = new Thread(() -> {
                while (true) {
                    String value;
                    while ((value = lines.poll()) != null) {
                        if (value.equals("END")) {
                            return;
                        }
                        pw.println(value);
                    }
                }
            });
            writeThread.start();

            final AtomicInteger counter = new AtomicInteger();

            cdr3AntigenMap.values().parallelStream().forEach(thisCdr3Info -> {
                        AminoAcidSequence thisCdr3 = thisCdr3Info.cdr3;
                        Cdr3Info otherCdr3Info;
                        NeighborhoodIterator<AminoAcidSequence, Cdr3Info> iter = stm
                                .getNeighborhoodIterator(thisCdr3, tsp);

                        Map<AminoAcidSequence, AlignmentInfo> alignmentVariants = new HashMap<>();

                        while ((otherCdr3Info = iter.next()) != null) {
                            Alignment<AminoAcidSequence> alignment = iter.getCurrentAlignment();

                            if (thisCdr3Info.nonDuplicateComparison(otherCdr3Info) &&
                                    alignment.getSequence1Range().length() == thisCdr3.size() &&
                                    !alignmentVariants.containsKey(otherCdr3Info.cdr3)) { // only one alignment
                                final AlignmentInfo alignmentInfo = new AlignmentInfo(
                                        thisCdr3Info.antigensOverlap(otherCdr3Info),
                                        alignment);

                                alignmentVariants.put(otherCdr3Info.cdr3, alignmentInfo);
                            }
                        }

                        for (Map.Entry<AminoAcidSequence, AlignmentInfo> alignmentEntry :
                                alignmentVariants.entrySet()) {

                            Mutations mutations = alignmentEntry.getValue().alignment.getAbsoluteMutations();

                            // Count true number of mismatches

                            int subst = 0, ins = 0, del = 0;

                            for (int k = 0; k < mutations.size(); k++) {
                                switch (mutations.getTypeByIndex(k)) {
                                    case Substitution:
                                        subst++;
                                        break;
                                    case Insertion:
                                        ins++;
                                        break;
                                    case Deletion:
                                        del++;
                                        break;
                                }
                            }

                            AminoAcidSequence otherCdr3 = alignmentEntry.getKey();
                            String edgeId = thisCdr3 + " (pp) " + otherCdr3,
                                    weight = Integer.toString(Math.max(otherCdr3.size(), thisCdr3.size()) -
                                            subst - 3 * (ins + del));

                            String line =
                                    thisCdr3 + "\t" + otherCdr3 + "\t" +
                                            (alignmentEntry.getValue().sameAntigen ? 1 : 0) + "\t" +
                                            subst + "\t" + ins + "\t" + del + "\t" + edgeId + "\t" + weight;

                            lines.add(line);
                        }

                        int count = counter.incrementAndGet();

                        if (count % 100 == 0) {
                            System.out.println("[" + (new Date()) + "] " +
                                    "Queried " + count + " of " + cdr3AntigenMap.size() + " cdr3 sequences.");
                        }
                    }
            );

            lines.add("END");

            writeThread.join();

            System.out.println("[" + (new Date()) + "] " +
                    "Done. Queried " + cdr3AntigenMap.size() + " cdr3 sequences.");
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static class Cdr3Info {
        final AminoAcidSequence cdr3;
        final Set<String> antigens = new HashSet<>();

        Cdr3Info(AminoAcidSequence cdr3) {
            this.cdr3 = cdr3;
        }

        void addAntigen(String antigen) {
            antigens.add(antigen);
        }

        boolean antigensOverlap(Cdr3Info other) {
            return overlaps(this.antigens, other.antigens);
        }

        boolean nonDuplicateComparison(Cdr3Info otherCdr3Info) {
            return cdr3.compareTo(otherCdr3Info.cdr3) > 0;
        }

        private static boolean overlaps(Set<String> set1, Set<String> set2) {
            if (set1.size() > set2.size()) {
                return overlaps(set2, set1);
            }
            for (String value : set1) {
                if (set2.contains(value))
                    return true;
            }
            return false;
        }
    }

    private static class AlignmentInfo {
        final boolean sameAntigen;
        final Alignment<AminoAcidSequence> alignment;

        AlignmentInfo(boolean sameAntigen, Alignment<AminoAcidSequence> alignment) {
            this.sameAntigen = sameAntigen;
            this.alignment = alignment;
        }
    }
}
