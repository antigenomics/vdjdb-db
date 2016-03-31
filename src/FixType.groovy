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


enum FixType {
    NoFixNeeded(false, true, 0),
    FixAdd(true, true, 2), FixTrim(true, true, 1),
    FixReplace(true, true, 3), FailedReplace(true, false, 5),
    FailedBadSegment(false, false, 4),
    FailedNoAlignment(true, false, 6)

    final boolean fixAttempted, good
    final int rank

    FixType(boolean fixAttempted, boolean good, int rank) {
        this.fixAttempted = fixAttempted
        this.good = good
        this.rank = rank
    }
}