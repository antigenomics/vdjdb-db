"""The embed extractor's three load-bearing transforms, and the contracts ``vdjdb-web`` depends on.

Any one of these failing silently blanks ``/overview``: the Scala side matches each image with a
regex that does not cross newlines, deletes nothing itself, and reads the fragment by a fixed
filename. So they are asserted here against a fixture rather than discovered in production.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path("summary/MakeEmbedableHtml.py")


def _module():
    spec = importlib.util.spec_from_file_location("make_embedable", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FIXTURE = """<html><head><title>x</title></head><body>
<p>before the marker, must not appear</p>
<p>!summary_embed_start!</p>
<div class="container">
<pre class="r"><code>ggplot(df)</code></pre>
<pre><code>## [1] &quot;printed output&quot;</code></pre>
<table>
<tr><th>a</th></tr>
</table>
<p><img role="img" src="data:image/png;base64,iVBORw0KGgo" width="1152" /></p>
</div>
<p>!summary_embed_end!</p>
<p>after the marker, must not appear</p>
</body></html>
"""


@pytest.fixture
def fragment():
    return "".join(_module().extract(FIXTURE.splitlines(keepends=True)))


def test_only_the_marked_region_is_extracted(fragment):
    assert "before the marker" not in fragment
    assert "after the marker" not in fragment
    assert "<html" not in fragment and "<head" not in fragment and "<body" not in fragment


def test_divs_are_dropped_and_r_source_is_skipped(fragment):
    # vdjdb-web injects the fragment into its own layout; a stray <div> breaks it.
    assert "<div" not in fragment
    assert "ggplot(df)" not in fragment


def test_printed_output_is_unwrapped_and_tables_get_semantic_ui_classes(fragment):
    assert "printed output" in fragment and "&quot;" not in fragment
    assert 'class="ui unstackable single line celled stripped compact small table"' in fragment


def test_images_are_responsive_and_the_dead_width_rewrite_is_gone(fragment):
    # #460. The old script rewrote width="1152" -> width="672", but pandoc emits no width at all,
    # so the rewrite never fired. A style on the tag cannot silently stop matching.
    assert 'style="max-width:100%;height:auto"' in fragment
    assert 'width="672"' not in fragment


def test_base64_payloads_stay_on_one_line(fragment):
    # The Scala regex does not cross newlines.
    payloads = [n for n in fragment.splitlines() if "data:image/png;base64," in n]
    assert payloads
    assert all(n.rstrip().endswith(("/>", ">", "</p>")) for n in payloads)
