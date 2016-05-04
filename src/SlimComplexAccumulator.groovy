/**
 * Created by mikesh on 5/3/16.
 */
class SlimComplexAccumulator {
    static final List<String> COMPLEX_SLIM_ANNOT_COLS = [
            "gene",
            "cdr3",
            "species",
            "antigen.epitope",
            "antigen.gene",
            "antigen.species"
    ],                        SUMMARY_COLS = ["complex.id",
                                              "v.segm", "d.segm", "j.segm",
                                              "mhc.a", "mhc.b", "mhc.class",
                                              "reference.id", "vdjdb.score"]

    final Set<String> complexId = new HashSet<>(),
                      vSegm = new HashSet<>(),
                      dSegm = new HashSet<>(),
                      jSegm = new HashSet<>(),
                      mhcA = new HashSet<>(), mhcB = new HashSet<>(),
                      mhcClass = new HashSet<>(),
                      referenceId = new HashSet<>()

    int vdjdbScore = 0

    SlimComplexAccumulator() {
    }

    void append(List<String> splitLine, Map<String, Integer> colIdMap) {
        complexId.addAll(splitLine[colIdMap["complex.id"]].split("[,;]").findAll { it.length() > 0 && it != "0" })
        vSegm.addAll(splitLine[colIdMap["v.segm"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        dSegm.addAll(splitLine[colIdMap["d.segm"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        jSegm.addAll(splitLine[colIdMap["j.segm"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcA.addAll(splitLine[colIdMap["mhc.a"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcB.addAll(splitLine[colIdMap["mhc.b"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcClass.addAll(splitLine[colIdMap["mhc.class"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        referenceId.addAll(splitLine[colIdMap["reference.id"]].split("[,;]").findAll { it.length() > 0 && it != "." })

        vdjdbScore = Math.max(vdjdbScore, splitLine[colIdMap["vdjdb.score"]].toInteger())
    }

    String getSummary() {
        [complexId.empty ? "0" : complexId.join(","),
         vSegm.join(","),
         dSegm.join(","),
         jSegm.join(","),
         mhcA.join(","),
         mhcB.join(","),
         mhcClass.join(","),
         referenceId.join(","),
         vdjdbScore].join("\t")
    }
}
