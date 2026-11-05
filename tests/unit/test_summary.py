"""The fragment `vdjdb-web` serves, produced by pandoc rather than by string surgery.

Any one of these failing silently blanks ``/overview``: the Scala side matches each image with a
regex that does not cross newlines, deletes nothing itself, and reads the fragment by a fixed
filename. They used to be properties of a Python script that line-scanned pandoc's output; they are
now properties of ``summary/embed.lua`` and ``summary/embed.tpl``, so the test drives the real
pipeline -- pandoc, the template and the filter -- over a fixture rather than a reimplementation.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
import re
import json

import pytest

from vdjdb.summary import render

FIXTURE = """---
title: fixture
---

before the marker, must not appear

!summary_embed_start!

#### A heading

``` r
ggplot(df) + geom_point()
```

```
## [1] "printed output"
```

<table>
<tr><th>a</th></tr>
</table>

<img src="fig.png" alt="" width="1152" />

!summary_embed_end!

after the marker, must not appear
"""

pandocmark = pytest.mark.skipif(shutil.which("pandoc") is None, reason="needs pandoc")


@pytest.mark.skipif(shutil.which("Rscript") is None or shutil.which("pandoc") is None,
                    reason="needs dashboard R/pandoc runtime")
def test_covid_plot_excludes_sars1_and_keeps_distinct_epitopes(tmp_path):
    source = Path("summary/vdjdb_summary.Rmd").read_text()
    block = next(b for b in re.findall(r"```\{r[^\n]*\}\n(.*?)```", source, re.S)
                 if "-> df.c" in b)
    table_block = next(b for b in re.findall(r"```\{r[^\n]*\}\n(.*?)```", source, re.S)
                       if "kable" in b and "df.c %>%" in b)
    script = tmp_path / "covid.R"
    script.write_text('''
library(dplyr); library(stringr); library(ggplot2); library(ggalluvial); library(knitr)
df <- data.frame(species="HomoSapiens", antigen.species=rep(c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-1"), each=30),
 antigen.gene="S", mhc.a="HLA-DPA1*01:03", mhc.b="HLA-DPB1*04:01", mhc.class="MHCII",
 antigen.epitope=rep(c("ALWAAAAAA", "ALWCCCCCC", "SARS1ONLY"), each=30),
 reference.id="PMID:1", gene="TRB")
df <- rbind(df, transform(df, gene="TRA"))
df.context <- transform(df, id=seq_len(nrow(df)))
''' + 'plot <- eval(parse(text=' + json.dumps(block) + '))\n' + '''
stopifnot(setequal(df.c$Epitope, c("ALWAAAAAA", "ALWCCCCCC")))
strata <- ggplot_build(plot)$data[[2]]
stopifnot(setequal(as.character(strata$stratum[strata$x == 3]), c("ALWAAAAAA", "ALWCCCCCC")))
''' + 'writeLines(as.character(eval(parse(text=' + json.dumps(table_block) + '))), ' +
        json.dumps(str(tmp_path / "table.html")) + ')\n')
    result = subprocess.run(["Rscript", str(script)], capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, result.stderr
    import html
    intermediate = tmp_path / "hla.knit.md"
    intermediate.write_text("!summary_embed_start!\n\n" +
                            (tmp_path / "table.html").read_text() +
                            "\n\n!summary_embed_end!\n")
    fragment = tmp_path / "hla.html"
    render.extract(intermediate, fragment)
    assert "HLA-DPA1*01/HLA-DPB1*04" in html.unescape(fragment.read_text())


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs dashboard R runtime")
def test_summary_keeps_observed_mhc_and_v_gene_pairs(tmp_path):
    import polars as pl
    from vdjdb.emit.legacy import SLIM_GROUP
    from vdjdb.summary.interactive import _panels

    rows = []
    for n in range(12):
        for a, b, v, ref in [("HLA-DPA1*01:03", "HLA-DPB1*04:01", "TRAV12-2*01", "PMID:1"),
                              ("HLA-DRA*01:01", "HLA-DRB1*04:05", "TRAV26-1*01", "PMID:2")]:
            rows.append(dict(zip([*SLIM_GROUP, "mhc.a", "mhc.b", "mhc.class", "v.segm", "reference.id"],
                                 ["TRA", "CAV" + "A" * (n + 1) + "F", "HomoSapiens", "GILGFVFTL", "M", "InfluenzaA", a, b, "MHCII", v, ref], strict=True)))
    long = pl.DataFrame(rows)
    slim = long.group_by(list(SLIM_GROUP)).agg(pl.exclude(SLIM_GROUP).unique().sort().str.join(","))
    long.write_csv(tmp_path / "vdjdb.txt", separator="\t")
    slim.write_csv(tmp_path / "vdjdb.slim.txt", separator="\t")
    source = Path("summary/vdjdb_summary.Rmd").read_text()
    chunks = re.findall(r"```\{r[^\n]*\}\n(.*?)```", source, re.S)
    load = next(b for b in chunks if "df = fread_generated" in b)
    mhc = next(b for b in chunks if "df.m =" in b)
    vh = next(b for b in chunks if "df.vhla =" in b).split("#df.vhla$v.segm.split")[0]
    script = tmp_path / "joint.R"
    script.write_text('params <- list(legacy=' + json.dumps(str(tmp_path)) + ')\n' +
                      load + '\n' + mhc + '\n' + vh + '''
stopifnot(nrow(df.m) == 2, all(df.m$Records == 12), all(df.m$Studies == 1))
stopifnot(nrow(df.vhla) == 2, all(df.vhla$records == 12))
stopifnot(all(df.vhla$v.segm.split[df.vhla$mhc.a.split == "HLA-DPA1*01"] == "TRAV12-2"))
stopifnot(all(df.vhla$v.segm.split[df.vhla$mhc.a.split == "HLA-DRA*01"] == "TRAV26-1"))
''')
    result = subprocess.run(["Rscript", str(script)], capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, result.stderr
    panels = _panels()
    contexts = panels.joint_contexts(tmp_path, slim)
    cells = panels.v_hla(contexts, min_records=1)
    assert contexts["id"].n_unique() == 12
    assert cells.height == 2
    assert cells["records"].to_list() == [12, 12]
    assert set(cells.select("mhc", "v").iter_rows()) == {
        ("DPA1*01 / DPB1*04", "TRAV12-2"), ("DRA*01 / DRB1*04", "TRAV26-1")}


@pytest.fixture
def fragment(tmp_path):
    """Built by `render.extract()` itself, not by a reimplementation of its pandoc call.

    The first version of this test assembled its own argv with absolute `--template` and
    `--lua-filter` paths. It passed while the real code path failed on a clean runner, because
    pandoc resolves a bare `--template NAME` from its DATA directory rather than the working
    directory. A test that rebuilds the command under test cannot catch a bug in the command.
    """
    src = tmp_path / "doc.knit.md"
    src.write_text(FIXTURE)
    # A 1x1 PNG, so `--embed-resources` has something real to inline.
    (tmp_path / "fig.png").write_bytes(bytes.fromhex(
        "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c4"
        "890000000a49444154789c6300010000050001" "0d0a2db4" "0000000049454e44ae426082"))
    out = tmp_path / "fragment.html"
    render.extract(intermediate=src, fragment=out)
    return out.read_text()


@pandocmark
def test_only_the_marked_region_is_emitted(fragment):
    assert "before the marker" not in fragment
    assert "after the marker" not in fragment
    assert "<html" not in fragment and "<head" not in fragment and "<body" not in fragment


@pandocmark
def test_r_source_is_dropped_and_no_divs_are_emitted(fragment):
    # vdjdb-web injects the fragment into its own layout; a stray <div> breaks it. The old script
    # deleted them line by line; this pass never passes --section-divs, so none exist to delete.
    assert "geom_point" not in fragment
    assert "<div" not in fragment


@pandocmark
def test_printed_output_is_unwrapped_and_tables_get_semantic_ui_classes(fragment):
    assert "printed output" in fragment
    assert "&quot;" not in fragment and "##" not in fragment
    assert 'class="ui unstackable single line celled stripped compact small table"' in fragment


@pandocmark
def test_images_are_responsive_and_inlined_on_one_line(fragment):
    # #460: knitr emits width="1152" and pandoc's resource embedding drops it, which is why the old
    # pixel rewrite never fired -- it ran on the stage where the attribute was already gone.
    assert 'style="max-width:100%;height:auto"' in fragment
    assert "data:image/png;base64," in fragment
    payloads = [n for n in fragment.splitlines() if "data:image/png;base64," in n]
    assert payloads and all(n.rstrip().endswith(("/>", ">", "</p>")) for n in payloads)


@pandocmark
def test_missing_markers_fail_loudly_rather_than_publishing_nothing(tmp_path):
    src = tmp_path / "doc.knit.md"
    src.write_text("# no markers here\n\njust prose.\n")
    proc = subprocess.run(
        ["pandoc", src.name, "--from", render.READER, "--to", "html4", "--standalone",
         "--template", str(render.TEMPLATE), "--lua-filter", str(render.FILTER),
         "-o", "out.html"],
        cwd=tmp_path, capture_output=True, text=True, check=False)
    assert proc.returncode != 0
    assert "no blocks between" in proc.stderr


def test_the_study_count_reads_every_reference_and_not_only_the_first():
    """A slim row merged from several publications comma-joins them in `reference.id`.

    Keeping field 1 reported 527 of the 636 distinct non-blank references in the slim table, missing
    109. Rendering one build both ways moved human TRA 348 -> 401, human TRB 431 -> 483, mouse TRA
    71 -> 84 and mouse TRB 88 -> 151. Four tables on the page carried the undercount.
    `vdjdb.assemble.evidence.support_counts` counts distinct `reference.id` on `records`, where the
    column holds one value, and this is what keeps the R side reading the same quantity.
    """
    rmd = (Path(__file__).resolve().parents[2] / "summary" / "vdjdb_summary.Rmd").read_text()
    assert 'str_split_fixed(reference.id, ",", n = Inf)[,1]' not in rmd
    assert rmd.count('length(unique(unlist(strsplit(reference.id, ","))))') == 5


def test_transparent_dashboard_keeps_white_reference_checks(tmp_path):
    """Alpha must not turn empty background into apparent black ink or an SSIM defect."""
    import base64
    import io

    Image = pytest.importorskip("PIL.Image")
    import importlib.util

    spec = importlib.util.spec_from_file_location("check_summary", Path("summary/check_summary.py"))
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)

    transparent = Image.new("RGBA", (20, 20), (0, 0, 0, 0))
    opaque = Image.new("RGB", (20, 20), "white")
    for image in (transparent, opaque):
        image.paste("black", (2, 2, 5, 5))
        image.paste("#e41a1c", (10, 10, 13, 13))

    def payload(image):
        buffer = io.BytesIO()
        image.save(buffer, format="PNG")
        return base64.b64encode(buffer.getvalue()).decode()

    encoded = payload(transparent)
    html = f'<img src="data:image/png;base64,{encoded}">'
    reference = tmp_path / "opaque.html"
    reference.write_text(f'<img src="data:image/png;base64,{payload(opaque)}">')
    assert checker._image(encoded).tobytes() == opaque.tobytes()
    assert checker.style(html) == []
    assert checker.perceptual(html, reference) == []
