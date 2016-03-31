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

class Table implements Iterable<Row> {
    final List<String> header
    final List<Row> rows
    final Map<String, Integer> indices = new HashMap<>()

    Table(String[] header, List<String[]> rows = []) {
        this.header = header.toList()
        this.rows = rows.collect { new Row(it) }

        header.eachWithIndex { String entry, int i ->
            indices[entry] = i
        }
    }

    void append(Table other) {
        other.rows.each { row ->
            def values = new String[header.size()]

            header.eachWithIndex { it, ind ->
                values[ind] = row[it]
            }

            rows << new Row(values)
        }
    }

    @Override
    Iterator<Row> iterator() {
        rows.iterator()
    }

    class Row {
        final List<String> values

        Row(String[] values) {
            this.values = new ArrayList<>(values.length)

            values.each { this.values << it }
        }

        String getAt(String columnName) {
            values[indices[columnName]]
        }

        String putAt(String columnName, String value) {
            values[indices[columnName]] = value
        }
    }
}
