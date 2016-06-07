/*
 * Copyright 2016 Mikhail Shugay
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


class VdjdbScoreFactory {
    final Map<String, Integer> scoreMap = new HashMap<>()

    VdjdbScoreFactory(Table masterTable) {
        masterTable.each { row ->
            def sign = getSignature(row)//, score = getScore(row)

            def score = 0
            if (row["meta.structure.id"].trim().length() > 0) {
                score += 7 // we have structure, any questions? :)
            } else {
                def verifyMethod = row["method.verification"].toLowerCase()

                if (verifyMethod.contains("targets")) {
                    // Verification with target cells
                    score += 3
                } else if (verifyMethod.contains("stain") || verifyMethod.contains("sort")) {
                    // Verification by cloning & re-staining
                    score += 2
                }
				
                def method = row["method.identification"].toLowerCase().trim()
                if (method.length() > 0) {
                    score += 1

                    def freqThreshold = method.contains("sort") ? 0.1 : 0.5

                    // If identified using tetramer use frequency to assign score
                    // Otherwise score is 0
                    def freq = row["method.frequency"].trim()

                    if (freq.length() > 0) {
                        def x = freq.split("/+").collect { it.toInteger() }

                        if (x[0] > 1 && (x[0] / (double) x[1]) >= freqThreshold) {
                            score += 1
                        }
                    }

                    if (row["method.singlecell"].trim().length() > 0) {
                        // Single-cell sequencing performed
                        score += (verifyMethod.contains("direct") ? 3 : 1) // direct verification by single-cell sorting of pMHC binding T-cells
                    }
                }
            }

            scoreMap[sign] = Math.max(scoreMap[sign] ?: 0, score)
        }
    }

    static final List<String> SIGNATURE_COLS = [
            "cdr3.alpha",
            "v.alpha",
            "j.alpha",
            "cdr3.beta",
            "v.beta",
            "j.beta",
            "species",
            "mhc.a",
            "mhc.b",
            "mhc.class",
            "antigen.epitope"
    ]

    int getScore(Table.Row row) {
        scoreMap[getSignature(row)] ?: 0
    }

    String getSignature(Table.Row row) {
        SIGNATURE_COLS.collect { row[it] }.join("\t")
    }
}
