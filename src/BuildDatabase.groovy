import groovy.json.JsonBuilder
import groovy.json.JsonSlurper

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
// Misc utils and classes, temporary files
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def err = { String message ->
    System.err.println(message)
    System.exit(1)
}

new File("../tmp/").mkdirs()

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
    (str == "" || str.length() >= 3) && // if CDR3 is specified, it should be at least 3AA long
            str =~ /^[ARNDCQEGHILKMFPSTWYVX]+$/
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

def antigenSpeciesDict = new HashMap<String, String>(),
    antigenGeneDict = new HashMap<String, String>()

new File("../patches/antigen_epitope_species_gene.dict").readLines().drop(1).forEach { 
    def splitLine = it.split("[\t ]+")
    antigenSpeciesDict[splitLine[0]] = splitLine[1]
    antigenGeneDict[splitLine[0]] = splitLine[2]
}

assert antigenSpeciesDict["NLVPMVATV"] == "CMV"
assert antigenGeneDict["NLVPMVATV"] == "pp65"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read, validate and concatenate chunks
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def chunksToBuild = args.findAll { !it.startsWith("-") } //args.length > 0 ? args[-1].split(",") as List<String> : []

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


    // patching
    table.each { row ->
        def epi = row["antigen.epitope"]

        if (antigenSpeciesDict[epi]) {
            row["antigen.species"] = antigenSpeciesDict[epi]
            row["antigen.gene"] = antigenGeneDict[epi]
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

def cdr3Fixer = new Cdr3Fixer("../res/segments.txt", "../res/segments.aaparts.txt")

println "Fixing CDR3 sequences (stage I)"

masterTable.addCol "cdr3fix.alpha"
masterTable.addCol "cdr3fix.beta"

masterTable.each { row ->
    ["alpha", "beta"].each { gene ->
        def cdr3 = row["cdr3.$gene"],
            vId = row["v.$gene"],
            jId = row["j.$gene"],
            species = row["species"]

        if (cdr3 != "") {
            vId = vId == "" ? cdr3Fixer.guessId(cdr3, species, gene, true) : vId
            row["v.$gene"] = vId
            jId = jId == "" ? cdr3Fixer.guessId(cdr3, species, gene, false) : jId
            row["j.$gene"] = jId

            def fixerResult = cdr3Fixer.fix(
                    cdr3,
                    vId,
                    jId,
                    species,
                    gene
            )

            row["cdr3.$gene"] = fixerResult.cdr3

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

// collapse complexes - TRA and TRB go on separate lines, complex id is used to group them to parent record

println "Writing flat table"

def METADATA_LINES = ["name\ttype\tvisible\tsearchable\tautocomplete\tdata.type\ttitle\tcomment",
                      "complex.id\ttxt\t0\t0\t0\tcomplex.id\tcomplex.id\tTCR alpha and beta chain records having the same complex identifier belong to the same T-cell clone.",
                      "gene\ttxt\t1\t1\t1\tfactor\tGene\tTCR chain: alpha or beta.",
                      "cdr3\tseq\t1\t1\t0\tcdr3\tCDR3\tTCR complementarity determining region 3 (CDR3) amino acid sequence.",
                      "v.segm\ttxt\t1\t1\t1\tfactor\tV\tTCR Variable segment allele.",
                      "j.segm\ttxt\t1\t1\t1\tfactor\tJ\tTCR Joining segment allele.",
                      "species\ttxt\t1\t1\t1\tfactor\tSpecies\tTCR parent species.",
                      "mhc.a\ttxt\t1\t1\t1\tfactor\tMHC A\tFirst MHC chain allele.",
                      "mhc.b\ttxt\t1\t1\t1\tfactor\tMHC B\tSecond MHC chain allele (defaults to Beta2Microglobulin for MHC class I).",
                      "mhc.class\ttxt\t1\t1\t1\tfactor\tMHC class\tMHC class (I or II).",
                      "antigen.epitope\tseq\t1\t1\t1\tpeptide\tEpitope\tAmino acid sequence of the epitope.",
                      "antigen.gene\ttxt\t1\t1\t1\tfactor\tEpitope gene\tRepresentative parent gene of the epitope.",
                      "antigen.species\ttxt\t1\t1\t1\tfactor\tEpitope species\tRepresentative parent species of the epitope.",
                      "reference.id\ttxt\t1\t1\t1\turl\tReference\tPubmed reference / URL / or submitter details in case unpublished.",
                      "method\ttxt\t1\t0\t0\tmethod.json\tMethod\tDetails on method used to assay TCR specificity.",
                      "meta\ttxt\t1\t0\t0\tmeta.json\tMeta\tVarious meta-information: cell subset, donor status, etc.",
                      "cdr3fix\ttxt\t1\t0\t0\tfixer.json\tCDR3fix\tDetails on CDR3 sequence fixing (if applied) and consistency between V, J and reported CDR3 sequence.",
                      "vdjdb.score\ttxt\t1\t1\t0\tuint\tScore\tVDJdb confidence score, the higher is the score the more confidence we have in the antigen specificity annotation of a given TCR clonotype/clone. Zero score indicates that there are insufficient method details to draw any conclusion.",
                      "web.method\ttxt\t0\t0\t1\t0\tfactor\tInternal",
                      "web.method.seq\ttxt\t0\t0\t1\t0\tfactor\tInternal",
                      "web.cdr3fix.nc\ttxt\t0\t0\t1\t0\tfactor\tInternal",
                      "web.cdr3fix.unmp\ttxt\t0\t0\t1\t0\tfactor\tInternal"],
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

/* Internal stuff for fast table filtering using VDJdb server */

def getWebMethod = { row ->
    def data = row["method.identification"].toLowerCase()

    if (data.contains("sort")) return "sort"
    else if (data.contains("culture") || data.contains("cloning") || data.contains("targets")) return "culture"
    return "other"
}

def getWebMethodSeq = { row ->
    if (row["method.singlecell"].trim().length() > 0)
        return "singlecell"

    def data = row["method.sequencing"].toLowerCase()

    if (data.contains("sanger")) return "sanger"
    else if (data.contains("-seq")) return "amplicon"
    return "other"
}

/* end */

def complexIdCounter = 0
new File("../database/vdjdb.txt").withPrintWriter { pw ->
    pw.println(METADATA_LINES[1..-1].collect { it.split("\t")[0] }.join("\t"))
    masterTable.each { row ->
        def complexAnnot = COMPLEX_ANNOT_COLS.collect { row[it] }.join("\t")

        def methodAnnot = new JsonBuilder(METHOD_COLUMNS.collectEntries {
            [(it.split("method.")[1]): row[it]]
        }).toString()

        def metaAnnotMap = META_COLUMNS.collectEntries {
            [(it.split("meta.")[1]): row[it]]
        }

        metaAnnotMap << [("samples.found"): scoreFactory.getSamplesDetected(row)]
        metaAnnotMap << [("studies.found"): scoreFactory.getStudiesDetected(row)]

        def metaAnnot = new JsonBuilder(metaAnnotMap).toString()

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
                            row["j.$it"],
                            complexAnnot,
                            methodAnnot,
                            metaAnnot,
                            row["cdr3fix.$it"],
                            row["vdjdb.score"],
                            getWebMethod(row), getWebMethodSeq(row), // Internal
                            ".", "." // Internal. Placeholders for CDR3 fixer results
                ].join("\t"))
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Re-align gene segments in the final table using nucleotide-on-amino acid aligner
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if (!args.contains("--no2fix")) {
    println "Fixing CDR3 sequences (stage II)"
    println "(it may take a while...)"

    def cmd = ["python", "AlignBestSegments.py", "../database/vdjdb_full.txt", "../database/vdjdb.txt", "../res/segments.txt"]
    def proc = cmd.execute()
    proc.waitForProcessOutput(System.out, System.err)

    if (proc.exitValue() != 0) {
        throw new RuntimeException("AlignBestSegments failed")
    }
}

def vdjdbLines = new File("../database/vdjdb.txt").readLines().collect { it.split("\t") as List }

def cdr3FixIdx = vdjdbLines[0].indexOf("cdr3fix"),
    webCdr3FixNcIdx = vdjdbLines[0].indexOf("web.cdr3fix.nc"),
    webCdr3FixUnmpIdx = vdjdbLines[0].indexOf("web.cdr3fix.unmp")

/* Internal stuff for fast table filtering using VDJdb server: process final cdr3fixer results */

new File("../database/vdjdb.txt").withPrintWriter { pw ->
    pw.println(vdjdbLines[0].join("\t"))
    def jsonSlurper = new JsonSlurper()
    vdjdbLines[1..-1].each { splitLine ->
        cdr3FixJson = jsonSlurper.parseText(splitLine[cdr3FixIdx])
        splitLine[webCdr3FixNcIdx] = cdr3FixJson.vCanonical && cdr3FixJson.jCanonical ? "no" : "yes"
        splitLine[webCdr3FixUnmpIdx] = cdr3FixJson.vEnd > -1 && cdr3FixJson.jStart > -1 ? "no" : "yes"
        pw.println(splitLine.join("\t"))
    }
}

/* end */

// Generate a slim version of database

println "Generating and writing slim database"

def SLIM_METADATA_LINES = [
        "name\ttype",
        "gene\ttxt",
        "cdr3\tseq",
        "species\ttxt",
        "antigen.epitope\tseq",
        "antigen.gene\ttxt",
        "antigen.species\ttxt",
        "complex.id\ttxt",
        "v.segm\ttxt",
        "j.segm\ttxt",
        "v.end\ttxt",
        "j.start\ttxt",
        "mhc.a\ttxt",
        "mhc.b\ttxt",
        "mhc.class\ttxt",
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
            colIdMap[it] = ind
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

