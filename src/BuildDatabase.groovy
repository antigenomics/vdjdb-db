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
    !str.startsWith("HLA") || str =~ /^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$/
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
            it.startsWith("PMID:") || it.startsWith("doi:") ||
                    it.startsWith("http://") || it.startsWith("https://") ||
                    it.toLowerCase().contains("unpublished")
        }
]

def correctors = ["antigen.species", "antigen.gene"].collectEntries {
    [(it): new File("../patches/$it").readLines().collect { it.split("\t") }]
}

def correct = { String text, String from, String to ->
    text.trim().equalsIgnoreCase(from) ? to : text
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read, validate and concatenate chunks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def chunksToBuild = args.length > 0 ? args[0].split(",") as List<String> : []

def chunkFiles = new File("../chunks/").listFiles().findAll { chunkFile ->
    def chunkName = chunkFile.name
    !chunkName.startsWith(".") && chunkName.endsWith(".txt") &&
            (chunksToBuild.empty || chunksToBuild.any { chunkName.toLowerCase().contains(it.toLowerCase()) })
}

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

    table.each { row ->
        correctors.each { Map.Entry<String, List<String[]>> replacer ->
            def val = row[replacer.key]
            replacer.value.each { fromto ->
                val = val.split(",").collect { correct(it, fromto[0], fromto[1]) }.join(",")
            }
            row[replacer.key] = val
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

            row["cdr3.$it"] = fixerResult.cdr3

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

def METADATA_LINES = ["name\ttype\tvisible\tsearchable\tdata.type\ttitle\tcomment",
                      "complex.id\ttxt\t0\t0\tcomplex.id\tcomplex.id\tTCR alpha and beta chain records having the same complex identifier belong to the same T-cell clone.",
                      "gene\ttxt\t1\t1\tfactor\tGene\tTCR chain: alpha or beta.",
                      "cdr3\tseq\t1\t1\tcdr3\tCDR3\tTCR complementarity determining region 3 (CDR3) amino acid sequence.",
                      "v.segm\ttxt\t1\t1\tfactor\tV\tTCR Variable segment identifier.",
                      "d.segm\ttxt\t1\t0\tfactor\tD\tTCR Diversity segment identifier (beta chain only).",
                      "j.segm\ttxt\t1\t1\tfactor\tJ\tTCR Joining segment identifier.",
                      "species\ttxt\t1\t1\tfactor\tSpecies\tParent species of a given TCR.",
                      "mhc.a\ttxt\t1\t1\tfactor\tMHC.A\tFirst MHC chain identifier.",
                      "mhc.b\ttxt\t1\t1\tfactor\tMHC.B\tSecond MHC chain identifier (set to beta2microglobulin for MHC class I).",
                      "mhc.class\ttxt\t1\t1\tfactor\tMHC.class\tMHC class (I or II).",
                      "antigen.epitope\tseq\t1\t1\tpeptide\tAntigen.Epitope\tAmino acid sequence of the antigen peptide.",
                      "antigen.gene\ttxt\t1\t1\tfactor\tAntigen.Gene\tParent gene of the antigen peptide.",
                      "antigen.species\ttxt\t1\t1\tfactor\tAntigen.Species\tParent species of the antigen peptide.",
                      "reference.id\ttxt\t1\t1\turl\tReference\tPubmed reference / URL / or submitter details in case unpublished.",
                      "method\ttxt\t1\t0\tmethod.json\tMethod\tDetails on method used to assay TCR specificity.",
                      "meta\ttxt\t1\t0\tmeta.json\tMeta\tVarious meta-information: cell subset, donor status, etc.",
                      "cdr3fix\ttxt\t1\t0\tfixer.json\tcdr3fix\tDetails on CDR3 sequence fixing (if applied) and consistency between V, J and reported CDR3 sequence.",
                      "vdjdb.score\ttxt\t1\t1\tuint\tscore\tVDJdb confidence score, the higher is the score the more confidence we have in the antigen specificity annotation of a given TCR clonotype/clone. 0 score indicates that there are insufficient method details to draw any conclusion."],
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
        if (row["cdr3.alpha"] == "" || row["cdr3.beta"] == "") {
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

// Generate a slim version of database

println "Generating and writing slim database"

def SLIM_METADATA_LINES = [
        "name\ttype",
        "complex.id\ttxt",
        "gene\ttxt",
        "cdr3\tseq",
        "v.segm\ttxt",
        "d.segm\ttxt",
        "j.segm\ttxt",
        "species\ttxt",
        "mhc.a\ttxt",
        "mhc.b\ttxt",
        "mhc.class\ttxt",
        "antigen.epitope\tseq",
        "antigen.gene\ttxt",
        "antigen.species\ttxt",
        "reference.id\ttxt",
        "vdjdb.score\ttxt"
]

new File("../database/vdjdb.slim.meta.txt").withPrintWriter { pw ->
    pw.println(SLIM_METADATA_LINES.join("\n"))
}

def firstLine = true

def colIdMap = [:]

def slimAccumulatorMap = new HashMap<String, SlimComplexAccumulator>()

def getSignature = { splitLine ->
    SlimComplexAccumulator.COMPLEX_SLIM_ANNOT_COLS.collect {
        splitLine[colIdMap[it]]
    }.join("\t")
}

new File("../database/vdjdb.txt").splitEachLine("\t") { splitLine ->
    if (firstLine) {
        def requiredIds = [SlimComplexAccumulator.COMPLEX_SLIM_ANNOT_COLS,
                           SlimComplexAccumulator.SUMMARY_COLS].flatten()

        splitLine.eachWithIndex { it, ind ->
//            if (requiredIds.contains(it)) {
                colIdMap[it] = ind
//            }
        }

        firstLine = false
    } else {
        def signature = getSignature(splitLine)

        def accumulator = slimAccumulatorMap[signature]

        if (accumulator == null) {
            slimAccumulatorMap[signature] = (accumulator = new SlimComplexAccumulator())
        }

        accumulator.append(splitLine, colIdMap)
    }
}

new File("../database/vdjdb.slim.txt").withPrintWriter { pw ->
    pw.println([SlimComplexAccumulator.COMPLEX_SLIM_ANNOT_COLS,
                SlimComplexAccumulator.SUMMARY_COLS].flatten().join("\t"))

    slimAccumulatorMap.each {
        pw.println(it.key + "\t" + it.value.summary)
    }
}

