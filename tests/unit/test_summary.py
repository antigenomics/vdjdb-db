"""The fragment `vdjdb-web` serves, produced by pandoc rather than by string surgery.

Any one of these failing silently blanks ``/overview``: the Scala side matches each image with a
regex that does not cross newlines, deletes nothing itself, and reads the fragment by a fixed
filename. They used to be properties of a Python script that line-scanned pandoc's output; they are
now properties of ``summary/embed.lua`` and ``summary/embed.html``, so the test drives the real
pipeline -- pandoc, the template and the filter -- over a fixture rather than a reimplementation.
"""
from __future__ import annotations

import shutil
import subprocess

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
