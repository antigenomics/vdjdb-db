VDJdb
=====

.. raw:: html

   <div class="proj-intro">
     <div>
       <p class="proj-intro__eyebrow">T-CELL RECEPTOR ANTIGEN SPECIFICITY</p>
       <p class="proj-intro__lead">VDJdb is a curated database of T-cell receptor sequences of known
       antigen specificity. <strong>A chunk is one publication</strong> and a chunk row is one
       record reporting both chains, so two matching rows in two chunks are independent
       reports, never duplicates. Every shipped file is a projection of two tidy tables, and every
       column, vocabulary and score rule on this site is rendered from the build&rsquo;s own field
       registry.</p>
       <p class="proj-intro__links">
         <a href="getting-started.html">Get the data</a>
         <span>&middot;</span>
         <a href="standards/chunk-format.html">Chunk format</a>
         <span>&middot;</span>
         <a href="standards/columns.html">Every column</a>
         <span>&middot;</span>
         <a href="submission.html">Submit</a>
       </p>
     </div>
   </div>

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html">
       <h3>Getting started</h3>
       <p>What is in a release, where to download it, and how to build one yourself.</p>
     </a>
     <a class="proj-card" href="standards/chunk-format.html">
       <h3>The chunk format</h3>
       <p>One file per publication. Every complex, method and meta column a submission may carry.</p>
     </a>
     <a class="proj-card" href="standards/columns.html">
       <h3>Column reference</h3>
       <p>Every shipped table, rendered from the field registry while the page builds.</p>
     </a>
     <a class="proj-card" href="standards/confidence-score.html">
       <h3>Confidence score</h3>
       <p>0&ndash;3, and what each level asserts about the evidence.</p>
     </a>
     <a class="proj-card" href="denoising.html">
       <h3>Motifs and denoising</h3>
       <p>What the motif stage is for, and the two-stage rule that picks its parameters.</p>
     </a>
     <a class="proj-card" href="clustering.html">
       <h3>Clustering</h3>
       <p>Six algorithms, 252 configurations, one harness, and the scorecard of the shipped one.</p>
     </a>
   </div>

.. toctree::
   :hidden:
   :caption: Getting started

   getting-started
   submission

.. toctree::
   :hidden:
   :caption: Data standards

   standards/chunk-format
   standards/columns
   standards/confidence-score
   standards/cdr3-fixing
   standards/identity
   standards/airr-mapping
   outputs

.. toctree::
   :hidden:
   :caption: Motif inference

   denoising
   clustering

.. toctree::
   :hidden:
   :caption: Builds

   builds
   dashboard
