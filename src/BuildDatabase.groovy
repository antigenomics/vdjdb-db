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
],  SIGNATURE_COLS = [
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

def err = { String message ->
    System.err.println(message)
    System.exit(1)
}

def checkHeader = { String[] header ->
    if (!header)
        err("Empty file")

    def colSet = new HashSet(header.toList())

    if (header.length != colSet.size())
        err("Duplicate columns found: $header")

    def missingColumns = [COMPLEX_COLUMNS, METHOD_COLUMNS, META_COLUMNS].flatten().findAll {
        !colSet.contains(it)
    }

    if (!missingColumns.empty)
        err("The following columns are missing: $missingColumns")
}

class Table {
    final String[] header
    final List<Row> rows
    final Map<String, Integer> indices = new HashMap<>()

    Table(String[] header, List<String[]> rows) {
        this.header = header
        this.rows = rows.collect { new Row(it) }

        header.eachWithIndex { String entry, int i ->
            indices[entry] = i
        }
    }

    class Row {
        final String[] values

        Row(String[] values) {
            this.values = new String[values.length]

            for (int i = 0; i < values.length; i++) {
                this.values[i] = values[i].trim()
            }
        }

        String getAt(String columnName) {
            values[indices[columnName]]
        }
    }
}


def isAASeqValid = { String str ->
    str =~ /^[ARNDCQEGHILKMFPSTWYV]+$/
}

def isMhcValid = { String str ->
    !str.startsWith("HLA") || str =~ /^HLA-[A-Z]+\*\d{2}(:\d{2})?$/
}

def speciesList = ["homosapiens", "musmusculus",
                   "rattusnorvegicus", "macacamulatta"]

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

def chunkFiles = new File("../chunks_new/").listFiles().toList()

if (chunkFiles.empty)
    err("No database chunks to process")

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

            if (rowErrorMessages) {
                chunkErrorMessages[rowSignature] = rowErrorMessages.toString()
            }
        }
    }

    if (!chunkErrorMessages.empty) {
        println([SIGNATURE_COLS, "error.message"].flatten().join(" | "))
        chunkErrorMessages.each {
            println(it.key + " | " + it.value)
        }
        err("There were errors processing $chunkFile.name")
    }
}