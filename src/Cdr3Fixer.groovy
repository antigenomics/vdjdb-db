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

class Cdr3Fixer {
    final Map<String, Map<String, String>> segmentsByIdBySpecies = new HashMap<>()
    final int maxReplaceSize, minHitSize

    static Map<String, String> nomenclatureConversions = new File("../patches/nomenclature.conversions")
        .readLines()
        .findAll { !it.startsWith("#") }
        .collectEntries { def splitLine = it.split("[\t ]+"); [(splitLine[0]): splitLine[1].split(",")[0]] }

    Cdr3Fixer(String segmentsFileName, int maxReplaceSize = 1, int minHitSize = 2) {
        this.maxReplaceSize = maxReplaceSize
        this.minHitSize = minHitSize

        new File(segmentsFileName).splitEachLine('\t') { splitLine ->
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

        if (species.toLowerCase() == "homosapiens") {
          def conversion = nomenclatureConversions[id]
          if (conversion) {
            id = conversion
          }
        }

        if (segmentsById == null)
            return null

        // Doing an extensive search
        // with all possible combinations

        [id, Util.simplifySegmentName(id)].flatten().collect { String it ->
            [it,
             "$it*01".toString(),
             it.contains("-") ? [] : (1..100).collect { int i -> "$it-$i*01".toString() }]
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

    OneSideFixerResult fix(String cdr3, String id, String species, boolean fivePrime) {
        String closestId = getClosestId(species, id)
        String segmentSeq = getSegmentSeq(species, closestId)

        if (!segmentSeq) {
            return new OneSideFixerResult(fivePrime ? cdr3 : cdr3.reverse(), closestId, FixType.FailedBadSegment)
        }

        if (!fivePrime) {
            cdr3 = cdr3.reverse()
            segmentSeq = segmentSeq.reverse()
        }

        def scanner = new KmerScanner(segmentSeq, minHitSize)

        def hit = scanner.scan(cdr3)

        if (hit) {
            if (hit.startInSegment == 0) {
                if (hit.startInCdr3 == 0) {
                    return new OneSideFixerResult(cdr3, closestId, FixType.NoFixNeeded, hit.matchSize)
                } else {
                    return new OneSideFixerResult(cdr3.substring(hit.startInCdr3), closestId, FixType.FixTrim, hit.matchSize)
                }
            } else {
                if (hit.startInCdr3 == 0) {
                    return new OneSideFixerResult(segmentSeq.substring(0, hit.startInSegment) + cdr3, closestId, FixType.FixAdd, hit.matchSize)
                } else if (hit.startInCdr3 <= maxReplaceSize) {
                    return new OneSideFixerResult(segmentSeq.substring(0, hit.startInSegment) + cdr3.substring(hit.startInCdr3), closestId, FixType.FixReplace, hit.matchSize)
                } else {
                    return new OneSideFixerResult(cdr3, closestId, FixType.FailedReplace)
                }
            }
        } else {
            return new OneSideFixerResult(cdr3, closestId, FixType.FailedNoAlignment)
        }
    }

    FixerResult fix(String cdr3, String vId, String jId, String species) {
        def vResult = vId.split(",").collect {
            fix(cdr3, it, species, true)
        }.min { it.fixType.rank }

        def jResult = jId.split(",").collect {
            fix(vResult.cdr3, it, species, false)
        }.min { it.fixType.rank }

        def newCdr3 = jResult.cdr3.reverse()

        //
        //       ---- 4
        // ---- 4
        // 0123456789
        // CASSLPKLFF
        //
        // 10 - 4 = 6

        new FixerResult(newCdr3, cdr3, newCdr3 != cdr3,
                vResult.x,
                jResult.x < 0 ? jResult.x : newCdr3.length() - jResult.x,
                vResult.segmentId, vResult.fixType,
                jResult.segmentId, jResult.fixType)
    }
}