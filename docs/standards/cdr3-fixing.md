# CDR3 fixing


At this stage, a series of checks is performed for CDR3 sequence and reported V/J segments:

* In case of *canonical* (starting with conserved ``C`` and ending with ``F/W``) CDR3 sequences: checks if 5' and 3' germline parts match corresponding V/J segment sequences.
* In case of truncated CDR3 sequences: adds conserved ``C/F/W`` residues. Can add more missing residues in case a relatively large contiguous V/J germline match is present.
* In case excessive germline part is reported (e.g. ``FGXG`` instead of simply ``F`` at CDR3 3' part), excessive residues are removed.
* Can correct mismatches in V/J germline regions in case a reliable non-contiguous V/J match is found.

The main reason behind that is that current immune repertoire sequencing (RepSeq) data processing software reports *canonical* clonotype sequences, high number antigen-specific TCR sequences present in literature are reported inconsistently. The latter greatly complicates annotation of RepSeq data using known antigen-specific TCR sequences.

In case of good V/J germline matching and errors in CDR3 sequence, the final CDR3 sequence in the database is replaced by its fixed version. The following report of CDR3 fixer is placed under ``cdr3fix.alpha`` and ``cdr3fix.beta`` columns, e.g.

```json
{
	"fixNeeded":true,
	"good":false,
	"cdr3":"CASSQDVGTGGVFALYF",
	"cdr3_old":"CASSQDVGTGGVFALY",
	"jFixType":"FixAdd",
	"jId":"TRBJ1-6*01",
	"jCanonical":true,
	"jStart":14,
	"vFixType":"FailedBadSegment",
	"vId":null,
	"vCanonical":true,
	"vEnd":-1
	}
```

and

```json
{
	"fixNeeded":true,
	"good":true,
	"cdr3":"CASSLSRGGNQPQYF",
	"cdr3_old":"CASSLSRGGNQPQY",
	"jFixType":"FixAdd",
	"jId":"TRBJ1-5*01",
	"jCanonical":true,
	"jStart":9,
	"vFixType":"NoFixNeeded",
	"vId":"TRBV14*01",
	"vCanonical":true,
	"vEnd":4
}
```

Field descriptions:

field | description
------|-------------
``fixNeeded`` | ``true`` if corrected CDR3 sequence differs from the original one, ``false`` otherwise
``good`` | ``true`` if the fix can be applied, ``false`` if the fix cannot be applied due to bad V/J entry or no V/J matching
``cdr3`` | Fixed CDR3 sequence
``cdr3_old`` | Original CDR3 sequence
``jFixType`` | Type of fix applied to CDR3 J germline part
``jCanonical`` | ``true`` if CDR3 ends with ``F`` or ``W``, ``false`` otherwise
``jId``  | J segment identifier
``jStart``  | A 0-based index of first CDR3 amino acid that belongs to J segment
``vFixType`` | Type of fix applied to CDR3 V germline part
``vCanonical`` | ``true`` if CDR3 starts with ``C``, ``false`` otherwise
``vId`` | V segment identifier
``vEnd``  | A 0-based index of the last CDR3 amino acid of V segment plus one

> **Note:**

> Possible V and J fix types: ``NoFixNeeded``, ``FixAdd``, ``FixReplace``, ``FixTrim``, ``FailedReplace`` (too many mismatches), ``FailedBadSegment`` (bad segment entry), ``FailedNoAlignment`` (no alignment at all)
