/*
 * Copyright 2015 Mikhail Shugay
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

class Cdr3Fixer {
    final Map<String, Map<String, String>> segmentsByIdBySpecies = new HashMap<>()
    final int maxReplaceSize, minHitSize

    Cdr3Fixer(int maxReplaceSize = 1, int minHitSize = 2) {
        this.maxReplaceSize = maxReplaceSize
        this.minHitSize = minHitSize

        new File("segments.txt").splitEachLine('\t') { splitLine ->
            if (!splitLine[0].startsWith("#")) {
                def species = splitLine[0].toLowerCase(),
                    type = splitLine[2].toLowerCase(), id = splitLine[3],
                    refPoint = splitLine[4].toInteger(), seq = splitLine[5]

                def segmentsById = segmentsByIdBySpecies[species]

                if (segmentsById == null)
                    segmentsByIdBySpecies.put(species, segmentsById = new HashMap<String, String>())

                if (type.startsWith("v") || type.startsWith("j")) {
                    boolean isJSegment = type.startsWith("j")
                    seq = isJSegment ? seq.substring(0, refPoint + 4) : seq.substring(refPoint - 3)

                    segmentsById.put(id, Util.translateLinear(seq, isJSegment))
                }
            }
        }
    }

    String getClosestId(String species, String id) {
        def segmentsById = segmentsByIdBySpecies[species.toLowerCase()]

        if (segmentsById == null)
            return null

        // Doing an extensive search
        // with all possible combinations

        [id, Util.simplifySegmentName(id)].collect {
            [id, "$id*01".toString(), (1..9).collect { "$id-$it*01".toString() }]
        }.flatten().find {
            segmentsById.containsKey(it)
        }
    }

    String getSegmentSeq(String species, String id) {
        if (id == null)
            return null

        def segmentsById = segmentsByIdBySpecies[species.toLowerCase()]

        if (segmentsById == null)
            return null

        segmentsById[id]
    }

    OneSideFixerResult fix(String cdr3, String segmentSeq) {
        def scanner = new KmerScanner(segmentSeq, minHitSize)

        def hit = scanner.scan(cdr3)

        if (hit) {
            if (hit.startInSegment == 0) {
                if (hit.startInCdr3 == 0) {
                    return new OneSideFixerResult(cdr3, FixType.NoFixNeeded)
                } else {
                    return new OneSideFixerResult(cdr3.substring(hit.startInCdr3), FixType.FixTrim)
                }
            } else {
                if (hit.startInCdr3 == 0) {
                    return new OneSideFixerResult(segmentSeq.substring(0, hit.startInSegment) + cdr3, FixType.FixAdd)
                } else if (hit.startInCdr3 <= maxReplaceSize) {
                    return new OneSideFixerResult(segmentSeq.substring(0, hit.startInSegment) + cdr3.substring(hit.startInCdr3), FixType.FixReplace)
                } else {
                    return new OneSideFixerResult(cdr3, FixType.FailedReplace)
                }
            }
        } else {
            return new OneSideFixerResult(cdr3, FixType.FailedNoAlignment)
        }
    }

    FixerResult fix(String cdr3, String vId, String jId, String species) {
        cdr3 = Util.replaceNonAa(cdr3)

        vId = getClosestId(species, vId)
        jId = getClosestId(species, jId)

        String vSeq = getSegmentSeq(species, vId),
               jSeq = getSegmentSeq(species, jId)

        def vResult = vSeq ? fix(cdr3, vSeq) :
                new OneSideFixerResult(cdr3, FixType.FailedBadSegment),
            jResult = jSeq ? fix(vResult.cdr3.reverse(), jSeq.reverse()) :
                    new OneSideFixerResult(cdr3.reverse(), FixType.FailedBadSegment)

        new FixerResult(jResult.cdr3.reverse(), vId, vResult.fixType, jId, jResult.fixType)
    }
}
