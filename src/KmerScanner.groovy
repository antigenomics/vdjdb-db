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

class KmerScanner {
    final int minHitSize
    final Map<String, Integer> kmers = new HashMap<>();

    KmerScanner(String seq, int minHitSize = 2) {
        this.minHitSize = minHitSize

        for (int i = minHitSize; i <= seq.length(); i++) {
            for (int j = 0; j <= seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)
                kmers.put(kmer, j)
            }
        }
    }

    SearchResult scan(String seq) {
        def bestHit = (minHitSize..<seq.length()).collect { int i ->
            (0..<(seq.length() - i)).collect { int j ->
                String kmer = seq.substring(j, j + i)

                def hit = kmers[kmer]

                if (hit != null) {
                    new SearchResult(hit, j, kmer.length())
                } else {
                    new SearchResult(-1, j, -1)
                }
            }
        }.flatten().max { it.matchSize }

        bestHit.matchSize > 0 ? bestHit : null

        /*
        // iterate from largest window to smallest one
        for (int i = seq.length(); i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < seq.length() - i; j++) {
                String kmer = seq.substring(j, j + i)

                def hit = kmers[kmer]

                if (hit != null) {
                    return new SearchResult(hit,
                            j,
                            kmer.length()
                    )
                }
            }
        }

        null*/
    }
}