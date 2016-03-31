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

def chunkFiles = new File("../chunks_new/").listFiles().toList()

if (chunkFiles.empty)
    err("No database chunks to process")

def masterTable = new Table(ALL_COLS as String[])

chunkFiles.each { chunkFile ->
    def table = readChunk(chunkFile)

    def chunkErrorMessages = [:]
    def rowSignatures = new HashSet<String>()

    table.rows.each { row ->
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

            if (row["cdr3.alpha"].length() == 0 && row["cdr3.beta"].length()) {
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

def cdr3Fixer = new Cdr3Fixer()