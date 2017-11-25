/*
 * Copyright 2016-2017 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import groovy.json.JsonSlurper

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
