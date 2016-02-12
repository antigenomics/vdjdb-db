[![Build Status](https://travis-ci.org/antigenomics/vdjdb-db.svg?branch=master)](https://travis-ci.org/antigenomics/vdjdb-db)

# VDJDB: A curated database of T-cell receptor sequences of known antigen specificity

This repository hosts the submissions to database and scripts to check, fix and build the database itself.

To build database from submissions, run the ``src/build_db.py`` python script.

To query the database for your immune repertoire sample(s) use the [VDJdb](https://github.com/antigenomics/vdjdb) software.

A web-based GUI for the database [VDJdb-server](https://github.com/antigenomics/vdjdb-server) is currently under development.

## Submission guide

To submit previously published sequence follow the steps below:

* Create an issue(s) labeled as ``paper`` and named by the paper pubmed id, ``PMID:XXXXXXX``. Note that if paper is a meta-study, you can mark it as ``meta-paper`` and link issues for its references in a reply to this issue.

* Create new branch and add chunk(s) for corresponding papers named as ``PMID_XXXXXXX``. Don't forget to close/reference corresponding issues in the commit message.

* Create a pull request for the branch and check if it passes the CI build. If there are any issues, modify them by fixing/removing entries as necessary.

CI tests include table format checks and CDR3 sequence checks via [FixCDR3](https://github.com/antigenomics/fixcdr3).

To view the list of papers that were not yet processed go [here](https://github.com/antigenomics/vdjdb-db/labels/paper).

Don't forget to add corresponding record to ``citations.txt``. If the paper contains the data itself just add ``PMID:XXXXXXX(tab)PMID:XXXXXXX``, 
in case the paper also provides sequences from previously published studies add ``PMID:XXXXXXX(tab)PMID:XXXXXXX,PMID:YYYYYYY,PMID:ZZZZZZZ``.

## Database specification

Each database submission in ``chunks/`` folder should have the following header and columns:

| column header   | column description           | allowed/example values                                                      |
|-----------------|------------------------------|-----------------------------------------------------------------------------|
| record.id       | unique record id             | ``VDJDBRXXXXXXXX``, where ``X`` stands for a digit                          |
| complex.id      | complex identifier           | ``VDJDBCXXXXXXXX`` or ``VDJDBC00000000`` if record represents unpaired data |
| cdr3            | CDR3 sequence                | amino acid sequence                                                         |
| v.segm          | Variable segment id          | ``TRBV20-1``, ``TRBV20-1*01``, etc                                          |
| j.segm          | Joining segment id           | ``TRBJ2-1``, ``TRBJ2-1*01``, etc                                            |
| gene            | TCR gene                     | ``TRA`` or ``TRB``                                                          |
| species         | TCR species                  | ``HomoSapiens``, ``MusMusculus``, etc                                       |
| mhc.a           | first MHC chain              | ``HLA-A*0201``, etc                                                         |
| mhc.b           | second MHC chain             | ``B2M``, etc                                                                |
| mhc.type        | MHC class                    | ``MHCI`` or ``MHCII``                                                       |
| antigen         | antigen sequence             | amino acid sequence                                                         |
| antigen.gene    | parent gene of the antigen   | ``pp65``, ``BMLF1``, ``MBP``, etc	                                       |
| antigen.species | antigen species              | ``CMV``, ``EBV``, ``HomoSapiens``, etc                                      |
| method          | specificity inference method | ``tetramer``, ``pentamer``, ``restimulation``, etc                          |
| reference       | reference type               | e.g. ``pubmed``                                                             |
| reference.id    | reference id                 | e.g.``PMID:XXXXXXX``                                                        |

The resulting composite database fille will contain all those columns with an addition of ``comment`` column. The latter stores 
a summary of **all** additional columns in every chunk as a compact list of key-value pairs in JSON format. 
For example if chunk #1 contains columns

tissue       | cell type
-------------|-----------
 ``spleen``  | ``cd8``
 ``spleen``  | ``cd4``

and chunk #2 contains column

comment              |
---------------------|
 ``mutated antigen`` |
 `` ``               |

the resulting database will contain the following column in addition to the required ones:

comment                                      |
---------------------------------------------|
``{ "tissue":"spleen", "cell type":"cd8"  }``|
``{ "tissue":"spleen", "cell type":"cd4"  }``|
``{ "comment":"mutated antigen" }``          |
`` ``                                        |

Records with the same complex identifier represent TCR:peptide:MHC complexes, where both TCR alpha and beta chains are known.