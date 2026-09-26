--[[
  Produce the publishable fragment directly, instead of rendering a whole document and then
  deleting parts of it with string surgery.

      pandoc vdjdb_summary.knit.md --template summary/embed.html --lua-filter summary/embed.lua ...

  `MakeEmbedableHtml.py` line-scanned pandoc's output for `<div`, `<pre class="r">` and `<table>`.
  Every one of those is a guess about markup pandoc happens to emit today, and each guess silently
  blanks vdjdb-web's /overview when it stops matching. Here the same four transforms are stated on
  the document tree before it is written, so there is nothing to guess:

    1. keep only the blocks between the two markers;
    2. drop the R source blocks (`code_folding: hide` renders them, the fragment must not);
    3. unwrap printed R output -- `## "Last updated on ..."` becomes plain text, as it always has;
    4. add Semantic UI's classes to tables and a responsive style to images.

  The `<div>` stripping is gone entirely, because this pass does not pass `--section-divs` and so
  emits none. That removes the failure mode rather than defending against it.
--]]

local START, STOP = "!summary_embed_start!", "!summary_embed_end!"
local TABLE_CLASS = "ui unstackable single line celled stripped compact small table"
-- #460. knitr emits `width="672"`/`width="1152"` on these tags and pandoc's resource embedding
-- drops the attribute when it inlines the image, which is why the old pixel-width rewrite never
-- fired: it ran on the stage where the attribute no longer existed. A style cannot be dropped that
-- way, and it is responsive, which a fixed pixel width never was.
local IMG_STYLE = 'style="max-width:100%;height:auto" '

local function is_marker(block, text)
  return block.t == "Para" and pandoc.utils.stringify(block) == text
end

function Pandoc(doc)
  local out, inside = {}, false
  for _, block in ipairs(doc.blocks) do
    if is_marker(block, STOP) then
      inside = false
    elseif is_marker(block, START) then
      inside = true
    elseif inside then
      out[#out + 1] = block
    end
  end
  if #out == 0 then
    error(("no blocks between %s and %s -- the markers moved or were removed, and the fragment "
           .. "would have been published empty"):format(START, STOP))
  end
  return pandoc.Pandoc(out, doc.meta)
end

function CodeBlock(el)
  for _, class in ipairs(el.classes) do
    if class == "r" then
      return {}
    end
  end
  -- Printed R output: `## "Last updated on 26 September, 2026"` -> the sentence itself.
  return pandoc.CodeBlock((el.text:gsub("#", ""):gsub('"', "")), el.attr)
end

-- Replacement FUNCTIONS, not strings: `max-width:100%` contains a `%`, which Lua's `gsub` reads as
-- a capture reference in a replacement string and rejects. A function returns the text verbatim.
local function patch(text)
  text = text:gsub("<table>", function() return '<table class="' .. TABLE_CLASS .. '">' end)
  text = text:gsub("<img ", function() return "<img " .. IMG_STYLE end)
  return text
end

function RawBlock(el)
  if el.format ~= "html" then
    return el
  end
  return pandoc.RawBlock("html", patch(el.text))
end

-- knitr emits each figure as a bare `<img ...>` line, which pandoc parses as a RawInline inside a
-- Para -- NOT a RawBlock. Filtering only blocks left every image unstyled, and the structural
-- check caught it: 8 images, 0 with a style attribute.
function RawInline(el)
  if el.format ~= "html" then
    return el
  end
  return pandoc.RawInline("html", patch(el.text))
end
