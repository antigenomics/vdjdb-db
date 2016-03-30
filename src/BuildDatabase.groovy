/**
 * Created by mikesh on 3/30/16.
 */




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

    def missingColumns = [COMPLEX_COLUMNS, METHOD_COLUMNS, METHOD_COLUMNS].flatten().findAll {
        !colSet.contains(it)
    }

    if (!missingColumns.empty)
        err("The following columns are missing: $missingColumns")
}

def readChunk = { String fileName ->
    def header
    def rows = []

    new File(fileName).withReader { reader ->
        header = reader.readLine().split("\t")
        checkHeader(header)

        def line
        while (line = reader.readLine()) {
            def splitLine = line.split("\t")
            if (splitLine.length != header.length)
                err("Row $splitLine has wrong number of columns")
            rows << splitLine
        }
    }

    new Chunk(header, rows)
}

def chunks = new File("../chunks_new/").listFiles().toList()

if (chunks.empty)
    err("No database chunks to process")

