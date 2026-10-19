VDJdb documentation
===================

VDJdb records published observations of T-cell receptors recognising peptide-MHC complexes.

.. raw:: html

   <div class="proj-intro">
     <p class="proj-intro__lead">From a published observation to a traceable receptor:pMHC record.
     Sequences, restriction, experimental evidence and source reports stay linked.</p>
     <p class="proj-intro__links"><a href="https://vdjdb.com">Search VDJdb</a>
     <span aria-hidden="true">·</span>
     <a href="https://github.com/antigenomics/vdjdb-db/releases">Download a release</a></p>
   </div>

Explore and use the data
------------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html"><h3>First records</h3>
       <p>Download a release and read observations with Python.</p></a>
     <a class="proj-card" href="outputs.html"><h3>Choose a format</h3>
       <p>Compare primary tables, legacy files and AIRR exports.</p></a>
     <a class="proj-card" href="dashboard.html"><h3>Database overview</h3>
       <p>Read current counts, publication growth and repertoire summaries.</p></a>
   </div>

Interpret a record
------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="standards/columns.html"><h3>Columns and identifiers</h3>
       <p>Look up fields, controlled vocabularies and linked tables.</p></a>
     <a class="proj-card" href="standards/record-scoring-v3.html"><h3>Confidence and evidence</h3>
       <p>Understand v3 forecasts, missing evidence and validation.</p></a>
     <a class="proj-card" href="denoising.html"><h3>Sequence motifs</h3>
       <p>Follow TCRNET and TCREmp annotations and their quality criteria.</p></a>
   </div>

Contribute and maintain
-----------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="submission.html"><h3>Submit a publication</h3>
       <p>Prepare a chunk, preserve source observations and review QC.</p></a>
     <a class="proj-card" href="tenx-dextramer.html"><h3>Process a 10x assay</h3>
       <p>Turn contigs and reagent counts into reviewed specificity calls.</p></a>
     <a class="proj-card" href="builds.html"><h3>Build and release</h3>
       <p>Assemble the database, verify outputs and prepare a release.</p></a>
   </div>

Start with :doc:`standards/terminology` for the distinction between peptide provenance and
recognition of a peptide-MHC complex. The :doc:`standards/confidence-score` page defines the
legacy 0-3 score; :doc:`standards/record-scoring-v3` describes the separate v3 forecasts.

.. toctree::
   :hidden:
   :caption: Tutorials

   getting-started

.. toctree::
   :hidden:
   :caption: How-to guides

   submission
   tenx-dextramer
   builds
   build-integrity

.. toctree::
   :hidden:
   :caption: Reference

   standards/chunk-format
   standards/columns
   standards/confidence-score
   standards/record-scoring-v3
   standards/identity
   standards/airr-mapping
   standards/corpus
   outputs

.. toctree::
   :hidden:
   :caption: Explanation

   standards/terminology
   standards/cdr3-fixing
   denoising
   clustering
   dashboard
