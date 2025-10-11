ignore = True
ignore2 = False

with open("vdjdb_summary_embed.html", "w") as pw:
    with open("vdjdb_summary.html", "r") as f:
        for line in f:
            if "!summary_embed_end!" in line:
                ignore = True
                ignore2 = False

            if not ignore and "<pre class=\"r\"><code>" in line:
                ignore2 = True

            if not ignore and not ignore2 and not line.startswith("<div") and not line.startswith("</div"):
                if line.startswith("<pre><code>##"):
                    line = line.replace("#", "").replace("&quot;", "")

                line = line.replace("<table>",
                                    "<table class=\"ui unstackable single line celled stripped compact small table\">")
                line = line.replace('width="1152"', 'width="672"')

                pw.write(line)

            if "!summary_embed_start!" in line:
                ignore = False
                ignore2 = False

            if not ignore and ignore2 and "</code></pre>" in line:
                ignore2 = False
