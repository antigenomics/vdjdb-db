"""Extract the publishable fragment of the dashboard for ``vdjdb-web``'s ``/overview``.

    uv run python summary/MakeEmbedableHtml.py [rendered.html] [fragment.html]

`vdjdb-web`'s ``Database.getSummaryFile`` reads ``vdjdb_summary_embed.html`` by that exact name and
injects it into an Angular page, so three properties of the output are load-bearing and any one of
them silently blanks ``/overview``:

1. each image's base64 payload is a **single unbroken line** -- the Scala side matches it with a
   regex that does not cross newlines;
2. there are **no ``<div>`` wrappers**, which is why every line starting with ``<div`` is dropped;
3. tables carry Semantic UI's classes, which pandoc does not emit and this script injects.

``summary/check_summary.py`` asserts all three against the produced file.

**The width rewrite this script used to carry was dead.** It replaced ``width="1152"`` with
``width="672"``, but pandoc emits ``<img role="img" src="data:...">`` with no ``width`` attribute at
all -- measured on the shipped fragment, which contains zero ``width="..."`` occurrences across its
8 images. That is issue #460: the cartoons render at natural size because nothing ever constrained
them. The fix is a style on the tag, which cannot silently stop matching the way a literal width
could.
"""
from __future__ import annotations

import sys
from pathlib import Path

START, END = "!summary_embed_start!", "!summary_embed_end!"
TABLE_CLASS = "ui unstackable single line celled stripped compact small table"
#: #460: responsive instead of a hardcoded pixel width, so it cannot go stale as the figures change.
IMG_STYLE = 'style="max-width:100%;height:auto" '


def extract(lines: list[str]) -> list[str]:
    """The fragment between the two markers, with the three transforms applied."""
    out: list[str] = []
    inside = False
    skipping_code = False
    for line in lines:
        if END in line:
            inside = False
        if inside and '<pre class="r"><code>' in line:
            skipping_code = True
        if inside and not skipping_code and not line.startswith(("<div", "</div")):
            if line.startswith("<pre><code>##"):
                line = line.replace("#", "").replace("&quot;", "")
            line = line.replace("<table>", f'<table class="{TABLE_CLASS}">')
            line = line.replace("<img ", f"<img {IMG_STYLE}")
            out.append(line)
        if START in line:
            inside, skipping_code = True, False
        if inside and skipping_code and "</code></pre>" in line:
            skipping_code = False
    return out


def main(src: Path, dst: Path) -> int:
    fragment = extract(src.read_text().splitlines(keepends=True))
    dst.write_text("".join(fragment))
    return len(fragment)


if __name__ == "__main__":
    here = Path(__file__).parent
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "vdjdb_summary.html"
    target = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "vdjdb_summary_embed.html"
    print(f"{main(source, target)} lines -> {target}")
