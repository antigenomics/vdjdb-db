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

    final static List<String> MHC_MULTIMER = ["tetramer", "pentamer", "multimer"]

    VdjdbScoreFactory(Table masterTable) {
        masterTable.each { row ->
            def sign = getSignature(row), score = getScore(row)

            if (MHC_MULTIMER.any {
                row["method.identification"].toLowerCase().contains("$it-sort")
            }) {
                // If identified using tetramer use frequency to assign score
                // Otherwise score is 0
                def freq = row["method.frequency"]
                if (freq.length() > 0) {
                    def x = freq.split("/").collect { it.toInteger() }

                    if (x[0] > 1 && (x[0] / (double) x[1]) >= 0.1) {
                        score += 2
                    }
                }
            }

            def verifyMethod = row["method.verification"].toLowerCase()

            if (verifyMethod.contains("targets")) {
                // Verification with target cells
                score += 4
            } else if (verifyMethod.contains("stain")) {
                // Verification by cloning & re-staining
                score += 2
            }

            if (row["method.singlecell"] != "") {
                // Single-cell sequencing performed
                score += 1
            }

            scoreMap[sign] = score
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
            "antigen.epitope",
    ]

    int getScore(Table.Row row) {
        scoreMap[getSignature(row)] ?: 0
    }

    String getSignature(Table.Row row) {
        SIGNATURE_COLS.collect { row[it] }.join("\t")
    }
}
