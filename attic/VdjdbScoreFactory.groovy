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

class VdjdbScoreFactory {
    final Map<String, Integer> scoreMap = new HashMap<>()
    final Map<String, List<String>> publicationMap = new HashMap<>()

    VdjdbScoreFactory(Table masterTable) {
        masterTable.each { row ->
            try {
                def sign = getSignature(row)

                // Assign publications
                def pubList = publicationMap[sign]
                if (pubList == null) {
                    publicationMap.put(sign, pubList = new ArrayList<String>())
                }
                pubList << row["reference.id"]

                // Compute score
                def score
                if (row["meta.structure.id"].trim().length() > 0) {
                    score = 3 // we have structure, any questions? :)
                } else {
                    /*
                        - We are sure that we have correct TCR sequence (+0 or +2)
                          sanger - several cells sequenced (2+)
                          single-cell - sure automatically
                          amplicon-seq - frequency is higher than 0.01

                        - We have moderate confidence about T-cell specificity (+0 or +2)
                          sort-based method - frequency is higher than 0.1
                          culture-based method - frequency is higher than 0.5
                          limiting dilution/culture performed - frequency is higher than 0.5

                        - We have high confidence about T-cell specificity (+0 or +3..+5)
                          direct method - best
                          stimulation-based - medium
                          sort-based worst
                     */

                    def freqStr = row["method.frequency"].trim()
                    def freq = getFrequency(freqStr),
                        count = getNumberOfCells(freqStr)

                    assert freq <= 1.0f

                    // sequencing
                    def seqScore = 1

                    def singleCell = row["method.singlecell"].trim().toLowerCase()

                    if (singleCell.length() > 0 && singleCell != "no") {
                        seqScore = 3
                    } else {
                        switch (row["method.sequencing"].trim().toLowerCase()) {
                            case "sanger":
                                seqScore = count >= 2 ? 3 : 2
                                break
                            case "amplicon-seq":
                                seqScore = freq >= 0.01f ? 3 : 1
                                break
                        }
                    }

                    // Moderate confidence regarding specificity
                    def identifyMethod = row["method.identification"].toLowerCase().trim()

                    def specScore1 = 0

                    if (cultureBasedIdentification(identifyMethod)) {
                        specScore1 = freq >= 0.5f ? 1 : 0
                    } else {
                        if (sortBasedMethod(identifyMethod)) {
                            specScore1 = freq >= 0.05f ? 1 : 0
                        } else if (stimulationBasedMethod(identifyMethod)) {
                            specScore1 = freq >= 0.25f ? 1 : 0
                        }
                    }

                    // High confidence regarding specificity
                    def verifyMethod = row["method.verification"].toLowerCase()

                    def specScore2 = 0

                    if (directVerification(verifyMethod)) {
                        // A direct method to observe TCR specificity
                        specScore2 = 3
                    } else if (stimulationBasedMethod(verifyMethod)) {
                        // Verification by target lysis, etc
                        specScore2 = 2
                        seqScore = 3 // tcr was cloned
                    } else if (sortBasedMethod(verifyMethod)) {
                        // Verification by cloning & re-staining
                        specScore2 = 1
                        seqScore = 3 // tcr was cloned
                    }

                    score = Math.min(seqScore, specScore1 + specScore2)
                }

                scoreMap[sign] = Math.max(scoreMap[sign] ?: 0, score)
            } catch (Exception e) {
                throw new RuntimeException("Error: " + e + "\nin " +
                    row["reference.id"] + " " + row["cdr3.alpha"] + " " + row["cdr3.beta"])
            }
        }
    }

    static int getNumberOfCells(String freq) {
        if (freq.length() == 0) {
            return 0
        }
        freq.contains('/') ? freq.split("/+")[0].toInteger() : 0
    }

    static float getFrequency(String freq) {
        if (freq.length() == 0) {
            return 0
        }
        if (freq.contains('/')) {
            def x = freq.split("/+").collect { it.toFloat() }
            return x[0] / x[1]
        } else if (freq.endsWith('%') && freq.length() > 1) {
            return freq[0..-2].toFloat() / 100.0f
        } else {
            return freq.toFloat()
        }
    }

    static boolean sortBasedMethod(String method) {
        // Sorting, magnetic separation, etc
        method.contains("sort") || method.contains("beads") ||
                method.contains("separation") || method.contains("stain")
    }

    static boolean stimulationBasedMethod(String method) {
        // Indentification / verification by cloning with ag-loaded target cells, etc
        method.contains("targets")
    }

    static boolean cultureBasedIdentification(String method) {
        // Cultured prior to sequencing => frequency doesn't tell much
        method.contains("culture") || method.contains("cloning")
    }

    static boolean directVerification(String method) {
        // Verification by direct monitoring of pMHC binding
        method.contains("direct")
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

    int getSamplesDetected(Table.Row row) {
        (publicationMap[getSignature(row)] ?: []).size()
    }

    int getStudiesDetected(Table.Row row) {
        (publicationMap[getSignature(row)] ?: []).unique().size()
    }

    String getSignature(Table.Row row) {
        SIGNATURE_COLS.collect { row[it] }.join("\t")
    }
}
