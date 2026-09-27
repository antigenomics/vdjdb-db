# CDR3 fixing

At this stage the build checks each CDR3 sequence against the reported V and J segments:

* For *canonical* CDR3 sequences, those starting with a conserved ``C`` and ending with ``F/W``: checks whether the 5' and 3' germline parts match the corresponding V/J segment sequences.
* For truncated CDR3 sequences: adds the conserved ``C/F/W`` residues. Further missing residues are added when a relatively large contiguous V/J germline match is present.
* When an excessive germline part is reported (e.g. ``FGXG`` instead of just ``F`` at the CDR3 3' part), the excess residues are removed.
* Mismatches in the V/J germline regions are corrected when a reliable non-contiguous V/J match is found.

Repertoire sequencing (RepSeq) data processing software reports *canonical* clonotype sequences, while a high number of antigen-specific TCR sequences in the literature are reported inconsistently. Fixing brings both into the same form, so that RepSeq data can be annotated with database records.

When the V/J germline match is good and the CDR3 sequence contains errors, the database carries the fixed sequence in place of the original. The fixer's report is stored in the ``cdr3fix.alpha`` and ``cdr3fix.beta`` columns, e.g.

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
``fixNeeded`` | ``true`` if the corrected CDR3 sequence differs from the original one, ``false`` otherwise
``good`` | ``true`` if the fix can be applied, ``false`` if it cannot, due to a bad V/J entry or no V/J match
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

> The V and J fix types are ``NoFixNeeded``, ``FixAdd``, ``FixReplace``, ``FixTrim``, ``FailedReplace`` (too many mismatches), ``FailedBadSegment`` (bad segment entry) and ``FailedNoAlignment`` (no alignment at all).
