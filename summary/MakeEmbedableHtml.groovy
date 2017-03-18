def ignore = true, ignore2 = false
new File("vdjdb_summary_embed.html").withPrintWriter { pw ->
pw.println("<div class=\"col-lg-10 col-lg-offset-1\">")
new File("vdjdb_summary.html").eachLine { line ->
    if (line.contains("!summary_embed_end!")) {
        ignore = true
        ignore2 = false
    }
    if (!ignore && line.contains("<pre class=\"r\"><code>")) {
        ignore2 = true
    }
    if (!ignore && !ignore2 && !line.startsWith("<div") && !line.startsWith("</div")) {
        if (line.startsWith("<pre><code>##")) {
            line = line.replaceAll("#", "").replaceAll("&quot;", "")
        }
        pw.println(line.replaceAll("<table>", "<table class=\"table\">")
                       .replaceAll(/width\s*=\s*"\d+"/, "width=\"100%\""))
    }
    if (line.contains("!summary_embed_start!")) {
        ignore = false
        ignore2 = false
    }
    if (!ignore && ignore2 && line.contains("</code></pre>")) {
        ignore2 = false
    }
}
pw.println("</div>")
}
