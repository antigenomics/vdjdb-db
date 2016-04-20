import groovy.json.JsonBuilder

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Database specification
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def COMPLEX_COLUMNS = [
        "cdr3.alpha",
        "v.alpha",
        "j.alpha",
        "cdr3.beta",
        "v.beta",
        "d.beta",
        "j.beta",
        "species",
        "mhc.a",
        "mhc.b",
        "mhc.class",
        "antigen.epitope",
        "antigen.gene",
        "antigen.species",
        "reference.id"
],  METHOD_COLUMNS = [
        "method.identification",
        "method.frequency",
        "method.singlecell",
        "method.sequencing",
        "method.verification"
],  META_COLUMNS = [
        "meta.study.id",
        "meta.cell.subset",
        "meta.subject.cohort",
        "meta.subject.id",
        "meta.replica.id",
        "meta.clone.id",
        "meta.epitope.id",
        "meta.tissue",
        "meta.donor.MHC",
        "meta.donor.MHC.method",
        "meta.structure.id"
],  ALL_COLS = [
        COMPLEX_COLUMNS, METHOD_COLUMNS, META_COLUMNS
].flatten(),
    SIGNATURE_COLS = [
            "cdr3.alpha",
            "v.alpha",
            "j.alpha",
            "cdr3.beta",
            "v.beta",
            "d.beta",
            "j.beta",
            "species",
            "mhc.a",
            "mhc.b",
            "mhc.class",
            "antigen.epitope",
            "antigen.gene",
            "antigen.species",
            "reference.id",
            "meta.study.id",
            "meta.cell.subset",
            "meta.subject.cohort",
            "meta.subject.id",
            "meta.replica.id",
            "meta.clone.id",
            "meta.tissue"
    ]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Misc utils and classes
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def err = { String message ->
    System.err.println(message)
    System.exit(1)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Table utils
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def checkHeader = { String[] header ->
    if (!header)
        err("Empty file")

    def colSet = new HashSet(header.toList())

    if (header.length != colSet.size())
        err("Duplicate columns found: $header")

    def missingColumns = ALL_COLS.findAll {
        !colSet.contains(it)
    }

    if (!missingColumns.empty)
        err("The following columns are missing: $missingColumns")
}


def readChunk = { File chunkFile ->
    def header = null, rows = []

    chunkFile.withReader { reader ->
        header = reader.readLine().split("\t")
        checkHeader(header)

        def line
        while (line = reader.readLine()) {
            if (line.endsWith("\t")) {
                line += " "
            }

            def splitLine = line.split("\t")

            if (splitLine.length != header.length)
                err("Row $splitLine has wrong number of columns")

            rows << splitLine
        }
    }

    new Table(header, rows)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic validation utils
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def isAASeqValid = { String str ->
    str =~ /^[ARNDCQEGHILKMFPSTWYV]+$/
}

def isMhcValid = { String str ->
    !str.startsWith("HLA") || str =~ /^HLA-[A-Z]+\*\d{2}(:\d{2})?$/
}

def speciesList = ["homosapiens", "musmusculus", "rattusnorvegicus", "macacamulatta"]

def validators = [
        "cdr3.alpha"     : isAASeqValid,
        "v.alpha"        : { it.startsWith("TRAV") },
        "j.alpha"        : { it.startsWith("TRAJ") },
        "cdr3.beta"      : isAASeqValid,
        "v.beta"         : { it.startsWith("TRBV") },
        "d.beta"         : { it.startsWith("TRBD") },
        "j.beta"         : { it.startsWith("TRBJ") },
        "species"        : { speciesList.any { species -> it.toLowerCase() == species } },
        "mhc.a"          : isMhcValid,
        "mhc.b"          : isMhcValid,
        "mhc.class"      : { it == "MHCI" || it == "MHCII" },
        "antigen.epitope": isAASeqValid,
        "reference.id"   : {
            it.startsWith("PMID:") || it.startsWith("doi:") || it.toLowerCase().contains("unpublished")
        }
]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read, validate and concatenate chunks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def chunkFiles = new File("../chunks_new/").listFiles().findAll { !it.name.startsWith(".") && it.name.endsWith(".txt") }

if (chunkFiles.empty)
    err("No database chunks to process")

def masterTable = new Table(ALL_COLS as String[])

println "Assembling chunks"

chunkFiles.each { chunkFile ->
    println "Processing $chunkFile.name"
    def table = readChunk(chunkFile)

    def chunkErrorMessages = [:]
    def rowSignatures = new HashSet<String>()

    table.each { row ->
        def rowSignature = SIGNATURE_COLS.collect { row[it] }.join(" | ")

        if (rowSignatures.contains(rowSignature)) {
            chunkErrorMessages[rowSignature] = "duplicate"
        } else {
            def rowErrorMessages = []

            validators.each {
                def cellValue = row[it.key]
                if (cellValue.length() > 0 && !it.value(cellValue)) {
                    rowErrorMessages << "bad." + it.key
                }
            }

            if (row["cdr3.alpha"] == "" && row["cdr3.beta"] == "") {
                rowErrorMessages << "no.cdr3"
            }

            if (row["antigen.epitope"] == "") {
                rowErrorMessages << "no.antigen.seq"
            }

            if (row["mhc.a"] == "" || row["mhc.b"] == "") {
                rowErrorMessages << "no.mhc"
            }

            if (rowErrorMessages) {
                chunkErrorMessages[rowSignature] = rowErrorMessages.toString()
            }
        }
    }

    if (!chunkErrorMessages.isEmpty()) {
        println([SIGNATURE_COLS, "error.message"].flatten().join(" | "))
        chunkErrorMessages.each {
            println(it.key + " | " + it.value)
        }
        err("There were errors processing $chunkFile.name")
    }

    masterTable.append(table)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fix CDR3 sequences
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def cdr3Fixer = new Cdr3Fixer()

println "Fixing CDR3 sequences"

masterTable.addCol "cdr3fix.alpha"
masterTable.addCol "cdr3fix.beta"

masterTable.each { row ->
    ["alpha", "beta"].each {
        if (row["cdr3.$it"] != "") {
            def fixerResult = cdr3Fixer.fix(
                    row["cdr3.$it"],
                    row["v.$it"],
                    row["j.$it"],
                    row["species"]
            )

            if (fixerResult.good && fixerResult.fixNeeded) {
                row["cdr3.$it"] = fixerResult.cdr3
            }

            row.values << new JsonBuilder(fixerResult).toString()
        } else {
            row.values << ""
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute VDJdb score
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

println "Computing VDJdb scores"

def scoreFactory = new VdjdbScoreFactory(masterTable)

masterTable.addCol "vdjdb.score"

masterTable.each { row ->
    row.values << scoreFactory.getScore(row).toString()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Write final table
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

println "Writing table"

new File("../database/").mkdir()

new File("../database/vdjdb_full.txt").withPrintWriter { pw ->
    pw.println(masterTable.header.join("\t"))
    masterTable.each {
        pw.println(it.values.join("\t"))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Write collapsed table
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// collapse complexes

println "Writing flat table"

def METADATA_LINES = ["name\ttype\tvisible\tsearchable\tdata.type",
                      "complex.id\ttxt\t0\t0\tcomplex.id",
                      "gene\ttxt\t1\t1\tfactor",
                      "cdr3\tseq\t1\t1\tcdr3",
                      "v.segm\ttxt\t1\t1\tfactor",
                      "d.segm\ttxt\t1\t0\tfactor",
                      "j.segm\ttxt\t1\t1\tfactor",
                      "species\ttxt\t1\t1\tfactor",
                      "mhc.a\ttxt\t1\t1\tfactor",
                      "mhc.b\ttxt\t1\t1\tfactor",
                      "mhc.class\ttxt\t1\t1\tfactor",
                      "antigen.epitope\tseq\t1\t1\tpeptide",
                      "antigen.gene\ttxt\t1\t1\tfactor",
                      "antigen.species\ttxt\t1\t1\tfactor",
                      "reference.id\ttxt\t1\t1\turl",
                      "method\ttxt\t1\t0\tmethod.json",
                      "meta\ttxt\t1\t0\tmeta.json",
                      "cdr3fix\ttxt\t1\t0\tfixer.json",
                      "vdjdb.score\ttxt\t1\t1\tuint"],
    COMPLEX_ANNOT_COLS = [
            "species",
            "mhc.a",
            "mhc.b",
            "mhc.class",
            "antigen.epitope",
            "antigen.gene",
            "antigen.species",
            "reference.id"]

new File("../database/vdjdb.meta.txt").withPrintWriter { pw ->
    pw.println(METADATA_LINES.join("\n"))
}

def complexIdCounter = 0
new File("../database/vdjdb.txt").withPrintWriter { pw ->
    pw.println(METADATA_LINES[1..-1].collect { it.split("\t")[0] }.join("\t"))
    masterTable.each { row ->
        def complexAnnot = COMPLEX_ANNOT_COLS.collect { row[it] }.join("\t")

        def methodAnnot = new JsonBuilder(METHOD_COLUMNS.collectEntries {
            [(it.split("method.")[1]): row[it]]
        }).toString(),
            metaAnnot = new JsonBuilder(META_COLUMNS.collectEntries {
                [(it.split("meta.")[1]): row[it]]
            }).toString()

        def complexId
        if (row["cdr3.alpha"] == "") {
            complexId = 0
        } else {
            complexIdCounter += 1
            complexId = complexIdCounter
        }

        ["alpha", "beta"].each {
            if (row["cdr3.$it"] != "") {
                pw.println([complexId,
                            it == "alpha" ? "TRA" : "TRB",
                            row["cdr3.$it"],
                            row["v.$it"],
                            it == "alpha" ? "" : row["d.$it"],
                            row["j.$it"],
                            complexAnnot,
                            methodAnnot,
                            metaAnnot,
                            row["cdr3fix.$it"],
                            row["vdjdb.score"]
                ].join("\t"))
            }
        }
    }
}
