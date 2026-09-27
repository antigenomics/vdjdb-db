"""The reference corpus: documents, a vocabulary, postings, and two questions you can ask of them.

One document per reference, four token families over it, and tf-idf weights. Built to reproduce what
``vdjdb.com/refsearch/`` serves and to answer what it cannot: whether a receptor feature goes with an
antigen because of itself or because of something it travels with.

``docs/standards/corpus.md`` is the reference page; ``ROADMAP.md`` section 12, phase 17, is the design.
"""
