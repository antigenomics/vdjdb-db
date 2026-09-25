"""Turn enriched clonotypes into clusters: the Hamming-1 graph, its components, and a layout.

Stage two. :mod:`vdjdb.motifs.tcrnet` says *which* clonotypes are more connected than the background
explains; this says *to each other*. Within one ``(species, gene, epitope)`` an edge joins two
enriched clonotypes one substitution apart -- the same ball the enrichment was scored over -- and a
connected component of at least :data:`MIN_CLUSTER` members is a cluster.

The edges come from one batched ``seqtree.Index.search_batch`` over the enriched set against itself,
not from a Python double loop: at 1,702 clonotypes for GILGFVFTL alone a quadratic scan is 1.4M
comparisons per epitope (CLAUDE.md hard rule 8, rung 2).

**Cluster numbering is content-derived, not a counter.** Raw component numbers depend on vertex
insertion order and on the graph library's internals, so they move between releases and every
bookmarked ``vdjdb.com`` motif URL breaks. Components are ordered by size descending then by their
lexicographically smallest member, so a cluster whose membership is unchanged keeps its number
without anything being stored between builds (CLAUDE.md hard rule 9, ROADMAP section 8.8).
"""
from __future__ import annotations

import polars as pl

from ..config import SEED

#: Smallest component that becomes a cluster. Matches the shipped files, whose smallest ``csz`` is
#: 5 over 1,928 cids -- below that there is no motif to read off a logo.
MIN_CLUSTER = 5

#: Species initial + chain initial, the legacy ``cid`` prefix (``H.B.GILGFVFTL.1``).
_INITIAL = {"HomoSapiens": "H", "MusMusculus": "M", "RattusNorvegicus": "R", "MacacaMulatta": "Q"}


def _edges(seqs: list[str], scope: str) -> list[tuple[int, int]]:
    """Every pair of ``seqs`` within ``scope`` of each other, as index pairs ``i < j``.

    One index, one batched search. ``search_batch`` returns hits in query order, so the query's own
    position is a self-hit and is dropped by the ``i < j`` filter along with each pair's mirror.
    """
    import seqtree

    subs, ins, dels, total = (int(x) for x in scope.split(","))
    index = seqtree.Index.build(seqs)
    params = seqtree.SearchParams(subs, ins, dels, total)
    return [(i, h.ref_id)
            for i, hits in enumerate(index.search_batch(seqs, params))
            for h in hits if i < h.ref_id]


def _components(n: int, edges: list[tuple[int, int]]) -> list[int]:
    """Connected-component label per vertex, by union-find. ``n`` vertices, ``edges`` undirected."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return [find(i) for i in range(n)]


def _layout(n: int, edges: list[tuple[int, int]]) -> list[tuple[float, float]]:
    """Force-directed ``(x, y)`` per vertex, for the web's motif view.

    Display coordinates only -- nothing downstream reads them. Seeded from
    :data:`vdjdb.config.SEED`, because an unseeded layout is a different picture every build and
    would make the file's digest meaningless (CLAUDE.md hard rule 7).
    """
    import random

    import igraph

    random.seed(SEED)
    g = igraph.Graph(n=n, edges=edges)
    return [(float(x), float(y)) for x, y in g.layout_fruchterman_reingold()]


def _repr_allele(calls: pl.Series) -> str:
    """The modal allele of a cluster, ties broken lexicographically so the answer is not the
    frame's order. Empty when every member's call is empty."""
    counts = (calls.to_frame("c").filter(pl.col("c") != "")
              .group_by("c").len().sort(["len", "c"], descending=[True, False]))
    return counts["c"][0] if counts.height else ""


def clusters(enriched: pl.DataFrame, *, scope: str = "1,0,0,1",
             min_cluster: int = MIN_CLUSTER) -> pl.DataFrame:
    """Cluster every ``(species, gene, epitope)`` group of :func:`~vdjdb.motifs.tcrnet.enriched_clonotypes`.

    Returns one row per clustered clonotype with ``cid``, ``csz``, ``x``, ``y`` and the cluster's
    representative V and J alleles -- the shape ``cluster_members.txt`` needs, minus the record
    columns that :mod:`vdjdb.motifs.emit` joins back on.
    """
    out: list[pl.DataFrame] = []
    for (species, gene, epitope), grp in enriched.group_by(
            ["species", "gene", "antigen.epitope"], maintain_order=True):
        grp = grp.sort("junction_aa", "v_call", "j_call")
        seqs = grp["junction_aa"].to_list()
        edges = _edges(seqs, scope)
        labels = _components(len(seqs), edges)
        xy = _layout(len(seqs), edges)

        g = grp.with_columns(
            pl.Series("__label", labels),
            pl.Series("x", [p[0] for p in xy]), pl.Series("y", [p[1] for p in xy]),
        )
        sizes = g.group_by("__label").agg(pl.len().alias("csz"),
                                          pl.col("junction_aa").min().alias("__first"))
        # Size descending, then the smallest member: a number that follows the content, not the
        # order the vertices happened to arrive in.
        order = (sizes.filter(pl.col("csz") >= min_cluster)
                      .sort(["csz", "__first"], descending=[True, False])
                      .with_row_index("__n", offset=1))
        if not order.height:
            continue
        prefix = f"{_INITIAL.get(species, species[:1])}.{gene[-1]}.{epitope}"
        g = (g.join(order.select("__label", "csz", "__n"), on="__label", how="inner")
               .with_columns((pl.lit(prefix) + "." + pl.col("__n").cast(pl.Utf8)).alias("cid"))
               .drop("__label", "__n"))
        g = g.join(
            g.group_by("cid").agg(
                pl.col("v_call").map_batches(_repr_allele, returns_scalar=True).alias("v.segm.repr"),
                pl.col("j_call").map_batches(_repr_allele, returns_scalar=True).alias("j.segm.repr"),
            ), on="cid")
        out.append(g)

    if not out:
        return pl.DataFrame(schema={"species": pl.Utf8, "gene": pl.Utf8, "antigen.epitope": pl.Utf8,
                                    "junction_aa": pl.Utf8, "cid": pl.Utf8, "csz": pl.UInt32})
    return pl.concat(out, how="vertical").sort("species", "gene", "antigen.epitope", "cid",
                                               "junction_aa", "v_call", "j_call")
