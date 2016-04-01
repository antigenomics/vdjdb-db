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

class FixerResult {
    final String cdr3, vId, jId
    final boolean vCanonical, jCanonical, fixNeeded
    final FixType vFixType, jFixType

    FixerResult(String cdr3, boolean fixNeeded, String vId, FixType vFixType, String jId, FixType jFixType) {
        this.cdr3 = cdr3
        this.fixNeeded = fixNeeded
        this.vId = vId
        this.vFixType = vFixType
        this.jFixType = jFixType
        this.jId = jId
        this.vCanonical = cdr3.startsWith("C")
        this.jCanonical = cdr3.endsWith("F") || cdr3.endsWith("W")
    }

    boolean isGood() {
        vFixType.good && jFixType.good
    }
}
