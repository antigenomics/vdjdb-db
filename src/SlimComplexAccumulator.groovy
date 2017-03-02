import groovy.json.JsonSlurper

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
                                              "v.segm", "j.segm",
                                              "v.end", "j.start",
                                              "mhc.a", "mhc.b", "mhc.class",
                                              "reference.id", "vdjdb.score"]

    private static final JsonSlurper json = new JsonSlurper()

    final Set<String> complexId = new HashSet<>(),
                      vSegm = new HashSet<>(),
                      jSegm = new HashSet<>(),
                      mhcA = new HashSet<>(),
                      mhcB = new HashSet<>(),
                      mhcClass = new HashSet<>(),
                      referenceId = new HashSet<>()

    int vdjdbScore = 0, vEnd = -1, jStart = -1

    SlimComplexAccumulator() {
    }

    void append(List<String> splitLine, Map<String, Integer> colIdMap) {
        complexId.addAll(splitLine[colIdMap["complex.id"]].split("[,;]").findAll { it.length() > 0 && it != "0" })
        vSegm.addAll(splitLine[colIdMap["v.segm"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        jSegm.addAll(splitLine[colIdMap["j.segm"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcA.addAll(splitLine[colIdMap["mhc.a"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcB.addAll(splitLine[colIdMap["mhc.b"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        mhcClass.addAll(splitLine[colIdMap["mhc.class"]].split("[,;]").findAll { it.length() > 0 && it != "." })
        referenceId.addAll(splitLine[colIdMap["reference.id"]].split("[,;]").findAll { it.length() > 0 && it != "." })

        def fixData = json.parseText(splitLine[colIdMap["cdr3fix"]])

        int vEnd1 = (int) fixData["vEnd"].toInteger(),
            jStart1 = (int) fixData["jStart"].toInteger()

        if (vEnd < 0 || vEnd < vEnd1) {
            vEnd = vEnd1
        }

        if (jStart < 0 || jStart > jStart1) {
            jStart = jStart1
        }

        vdjdbScore = Math.max(vdjdbScore, splitLine[colIdMap["vdjdb.score"]].toInteger())
    }


    String getSummary() {
        [complexId.empty ? "0" : complexId.join(","),
         vSegm.join(","),
         jSegm.join(","),
         vEnd,
         jStart,
         mhcA.join(","),
         mhcB.join(","),
         mhcClass.join(","),
         referenceId.join(","),
         vdjdbScore].join("\t")
    }
}
