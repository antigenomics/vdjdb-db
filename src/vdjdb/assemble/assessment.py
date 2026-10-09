"""Reported peptides and allele-conditioned presentation predictions.

Predictions annotate observations; they never rewrite them or supply identity keys.
The reference panel is an explicit, checksum-verified input. No download or stored
calibration is used during assembly. Each supported MHC species/class is scored
once on its distinct peptide set and joined back to all reported provenance.
"""
from __future__ import annotations

import gc
import hashlib
import math
import os
import tomllib
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from copy import copy
from multiprocessing import get_context
from pathlib import Path

import polars as pl

from ..config import SEED, Paths
from ..schema import EPITOPE_ASSESSMENT_COLUMNS

PROVENANCE = ("antigen.epitope", "antigen.species", "antigen.gene", "species", "mhc.species", "mhc.class")
KEY = (*PROVENANCE, "mhc.a", "mhc.b")
_SPECIES = {"HomoSapiens": "human", "MusMusculus": "mouse"}
_CLASS = {"MHCI": "mhc1", "MHCII": "mhc2"}
_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
_THREAD_ENV = ("POLARS_MAX_THREADS", "RAYON_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "MKL_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "NUMEXPR_NUM_THREADS")


def specification(root: Path | None = None) -> dict:
    return tomllib.loads(((root or Paths.discover().root) / "rules" /
                          "epitope_assessment.toml").read_text())


def verify_reference(path: Path, spec: dict) -> None:
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != spec["sha256"]:
        raise ValueError("pMHC reference checksum differs from rules/epitope_assessment.toml")


def fetch_reference(out: Path) -> Path:
    """Fetch a pinned input, separately from a build. Existing bytes are verified."""
    from huggingface_hub import hf_hub_download

    spec = specification()
    path = Path(hf_hub_download(spec["repository"], spec["file"], repo_type="dataset",
                               revision=spec["revision"], local_dir=out))
    verify_reference(path, spec)
    return path


@contextmanager
def _environment(*, parallel: bool):
    # Disable dependency calibration persistence even in serial calls; restore caller settings.
    changes = {"MHCMATCH_CALIBRATION_CACHE": "off"}
    if parallel:
        changes.update(dict.fromkeys(_THREAD_ENV, "1"))
    previous = {key: os.environ.get(key) for key in changes}
    os.environ.update(changes)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _score_group(task: tuple) -> list[dict]:
    """Class-I scorer per species; allele-owned class-II scorers bound frame memoization."""
    import mhcmatch
    from mhcmatch.predict import KMER_LENS, band_for, build_scorer, tile
    from mhcmatch.pseudoseq import class2_key, resolve_allele
    from mhcmatch.store import binding_core

    path, spec, species, cls, observed = task[:5]
    shard, shards = task[5:] if len(task) > 5 else (0, 1)
    store = mhcmatch.Store.from_pmhc(path=path, species=species, classes=(cls,))
    panel = sorted(store.panel_alleles(cls))
    if not panel:
        return [{"antigen.epitope": p, "mhc.a": a, "mhc.b": b, "reported": True,
                 "assessment.status": "empty_panel"} for p, a, b in observed if shard == 0]
    # Resolve each distinct reported pair once. Do not impute an absent DP/DQ partner.
    aliases = {}
    for a, b in sorted({(a, b) for _, a, b in observed}):
        partner = "" if species == "mouse" and a == b else b
        name = class2_key(a, partner, impute_alpha=False) if cls == "mhc2" else a
        incomplete = cls == "mhc2" and any("DQ" in v or "DP" in v for v in (a, b)) and (not a or not b)
        hits = store.panel_alleles(cls, [name]) if name and not incomplete else []
        _resolved, exact = resolve_allele(name, cls) if name else (None, False)
        aliases[a, b] = (hits[0], "exact" if exact or name == hits[0] else "nearest") if hits else ("", "")
    rows = []
    selected = panel[len(panel) * shard // shards:len(panel) * (shard + 1) // shards]
    peptides = sorted({p for p, _, _ in observed})
    valid = []
    for peptide in peptides:
        reported = sorted({(a, b) for p, a, b in observed if p == peptide})
        if not set(peptide) <= _AA or len(peptide) < (8 if cls == "mhc1" else 9):
            rows.extend({"antigen.epitope": peptide, "mhc.a": a, "mhc.b": b,
                         "assessment.status": "unsupported_peptide", "reported": True,
                         "prediction.allele": aliases[a, b][0],
                         "allele.resolution": aliases[a, b][1]} for a, b in reported if shard == 0)
        else:
            valid.append(peptide)
    scored = {peptide: [] for peptide in valid}
    reference_store = store
    if valid and cls == "mhc1":
        model, cal, affinity = build_scorer(store, cls, background=spec["background"],
                                          footprint=spec["footprint"], seed=SEED)
    for allele in selected if valid else []:
        if cls == "mhc2":
            # Fresh public scorer ownership per allele bounds frame memoization over lengths.
            # Remove this workaround after antigenomics/mhcmatch#4 ships a bounded batch API.
            # Copy the unscored store: reference panels are read once and shared within this run.
            store = copy(reference_store)
            model, cal, affinity = build_scorer(store, cls, background=spec["background"],
                                              footprint=spec["footprint"], seed=SEED)
        for peptide in valid:
            candidates = list(tile(peptide, KMER_LENS[cls])) if cls == "mhc1" and len(peptide) > 11 \
                else [(peptide, 0)]
            best = None
            for candidate, offset in candidates:
                score = model.score(candidate, allele)
                if not math.isfinite(score):
                    continue
                rank = cal.percent_rank(allele, score, length=len(candidate) if cls == "mhc2" else None)
                if not math.isfinite(rank):
                    continue
                value = (rank, offset, candidate, score)
                if best is None or value[:3] < best[:3]:
                    best = value
            if best is None:
                continue
            rank, offset, candidate, score = best
            register = model.best_register(candidate, allele)[0] if cls == "mhc2" else None
            core, core_offset = binding_core(candidate, cls, register_start=register)
            face = store.decompose(candidate, cls, allele, register_start=register).tcr_facing
            core_face, _ = binding_core(face, cls, register_start=register)
            probability = cal.p_present(allele, score)
            scored[peptide].append((rank, allele, {
                "antigen.epitope": peptide, "prediction.allele": allele,
                "prediction.peptide": candidate, "prediction.offset": str(offset),
                "presentation.percent_rank": f"{rank:.3f}",
                "presentation.p_present": f"{probability:.4f}" if math.isfinite(probability) else "",
                "presentation.band": band_for(rank, cls), "core": core,
                "core.offset": str(core_offset) if core else "",
                "core.source": ("model" if cls == "mhc2" else "footprint") if core else "",
                "tcr.facing": face, "core.tcr.facing": core_face,
                "assessment.status": "scored",
                "__rank": rank,
            }))
        if cls == "mhc2":
            del model, cal, affinity, store
            gc.collect()
    for peptide in valid:
        reported = sorted({(a, b) for p, a, b in observed if p == peptide})
        predictions = scored[peptide]
        predictions.sort(key=lambda item: item[:2])
        found = set()
        for index, (_, allele, prediction) in enumerate(predictions):
            pairs = [(a, b) for a, b in reported if aliases[a, b][0] == allele]
            if not pairs and prediction["presentation.band"] == "non-binder" and index != 0:
                continue
            # Retain every reported pair even if it is a non-binder, plus all weak/strong
            # predictions and the panel's best allele even when no predicted presenter exists.
            for a, b in [*pairs, ("", "")]:
                rows.append({**prediction, "mhc.a": a, "mhc.b": b, "reported": (a, b) in pairs,
                             "prediction.best": index == 0,
                             "allele.resolution": aliases[a, b][1] if (a, b) in pairs else ""})
                found.add((a, b))
        rows.extend({"antigen.epitope": peptide, "mhc.a": a, "mhc.b": b,
                     "prediction.allele": aliases[a, b][0],
                     "allele.resolution": aliases[a, b][1],
                     "reported": True, "assessment.status": "not_scorable" if aliases[a, b][0]
                     else "allele_not_in_panel"}
                    for a, b in reported if (a, b) not in found and
                    (aliases[a, b][0] in selected or (not aliases[a, b][0] and shard == 0)))
    return rows


def build_assessment(records: pl.DataFrame, *, reference: Path | None = None,
                     jobs: int = 1) -> pl.DataFrame:
    """One row per provenance and reported pair or predicted panel allele.

    Optional measurements are strings with empty missing values, including in TSV.
    Numbers can be cast after selecting scored rows. Reported support is counted
    only on reported rows. Unsupported inputs remain visible, with explicit status.
    """
    if jobs < 1:
        raise ValueError("epitope jobs must be positive")
    schema = {c: pl.String for c in EPITOPE_ASSESSMENT_COLUMNS}
    schema.update({"reported": pl.Boolean, "prediction.best": pl.Boolean,
                   "records": pl.UInt32, "references": pl.UInt32})
    # MHC species follows the reported molecule, independently of receptor or antigen species.
    hla = pl.any_horizontal(pl.col(c).str.starts_with("HLA-") for c in ("mhc.a", "mhc.b"))
    mouse = pl.any_horizontal(pl.col(c).str.contains(r"^(H2-|H-2|I-[AE])")
                              for c in ("mhc.a", "mhc.b"))
    inputs = records.with_columns(
        pl.when(hla & ~mouse).then(pl.lit("HomoSapiens"))
        .when(mouse & ~hla).then(pl.lit("MusMusculus"))
        .otherwise(pl.lit("")).alias("mhc.species"))
    counts = (inputs.group_by(KEY).agg(pl.len().alias("records"),
                                      pl.col("reference.id").n_unique().alias("references"))
              .sort(KEY))
    if counts.is_empty():
        return pl.DataFrame(schema=schema)
    spec = specification()
    predictions = []
    if reference is not None:
        import mhcmatch

        verify_reference(reference, spec)
        if mhcmatch.__version__ != spec["mhcmatch_version"]:
            raise ValueError("mhcmatch version differs from rules/epitope_assessment.toml")
        unique = counts.select("mhc.species", "mhc.class", "antigen.epitope", "mhc.a", "mhc.b").unique()
        tasks = []
        groups = unique.partition_by("mhc.species", "mhc.class", as_dict=True)
        for (species, cls), frame in sorted(groups.items()):
            if species not in _SPECIES or cls not in _CLASS:
                continue
            observed = frame.select("antigen.epitope", "mhc.a", "mhc.b").sort(
                "antigen.epitope", "mhc.a", "mhc.b").rows()
            slices = jobs if cls == "MHCII" else 1
            tasks.extend((str(reference), spec, _SPECIES[species], _CLASS[cls],
                          observed, i, slices) for i in range(slices))
        with _environment(parallel=jobs > 1):
            if jobs > 1 and len(tasks) > 1:
                with ProcessPoolExecutor(max_workers=min(jobs, len(tasks)),
                                         mp_context=get_context("spawn")) as pool:
                    results = list(pool.map(_score_group, tasks))
            else:
                results = [_score_group(task) for task in tasks]
        for task, result in zip(tasks, results, strict=True):
            species = next(key for key, value in _SPECIES.items() if value == task[2])
            cls = next(key for key, value in _CLASS.items() if value == task[3])
            predictions.extend({**row, "mhc.species": species, "mhc.class": cls} for row in result)
    pred_schema = {c: t for c, t in schema.items() if c not in ("antigen.gene", "antigen.species",
                                                             "species", "records", "references")}
    # Expand predictions to every reported provenance, including competing parent-gene labels.
    defaults = {c: False if t == pl.Boolean else "" for c, t in pred_schema.items()}
    predicted = pl.DataFrame([{**defaults, "__rank": float("inf"), **row} for row in predictions],
                             schema={**pred_schema, "__rank": pl.Float64})
    group = ["antigen.epitope", "mhc.species", "mhc.class"]
    best = (predicted.filter(pl.col("assessment.status") == "scored")
            .sort("__rank", "prediction.allele").group_by(group)
            .agg(pl.col("prediction.allele").first().alias("__best")))
    predicted = (predicted.join(best, on=group, how="left")
                 .with_columns(((pl.col("assessment.status") == "scored") &
                                (pl.col("prediction.allele") == pl.col("__best")))
                               .fill_null(False).alias("prediction.best"))
                 .drop("__rank", "__best"))
    provenance = counts.select(PROVENANCE).unique()
    expanded = provenance.join(predicted, on=["antigen.epitope", "mhc.species", "mhc.class"], how="inner")
    # A prediction for an alias must not claim that every source reports that alias.
    expanded = (expanded.join(counts, on=KEY, how="left")
                .filter(~pl.col("reported") | pl.col("records").is_not_null()))
    reported_rows = expanded.filter(pl.col("reported"))
    inferred_rows = expanded.filter(
        ~pl.col("reported") & ((pl.col("presentation.band") != "non-binder")
                               | pl.col("prediction.best"))).join(
        reported_rows.filter(pl.col("prediction.allele") != "").select(
            *PROVENANCE, "prediction.allele").unique(),
        on=[*PROVENANCE, "prediction.allele"], how="anti")
    expanded = pl.concat([reported_rows, inferred_rows])
    covered = reported_rows.select(KEY)
    unscored = (counts.join(covered, on=KEY, how="anti")
                .with_columns(pl.lit(True).alias("reported"),
                              pl.lit("reference_not_supplied" if reference is None else
                                     "unsupported_mhc_species_or_class").alias("assessment.status")))
    result = pl.concat([expanded, unscored], how="diagonal")
    result = result.with_columns(*(pl.col(c).fill_null(False if t == pl.Boolean else
                                                      0 if t == pl.UInt32 else "").cast(t)
                                   for c, t in schema.items()))
    return (result.with_columns(pl.lit(spec["background"] if reference else "")
                                .alias("prediction.background"),
                                pl.lit(spec["footprint"] if reference else "").alias("prediction.footprint"),
                                pl.lit(spec["mhcmatch_version"] if reference else "")
                                .alias("mhcmatch.version"),
                                pl.lit(spec["sha256"] if reference else "").alias("reference.sha256"),
                                pl.lit(str(SEED) if reference else "").alias("prediction.seed"))
            .select(EPITOPE_ASSESSMENT_COLUMNS).sort(*KEY, "prediction.allele"))
