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

-- Set with `-M asset-prefix=...` to reference the figures instead of inlining them. Empty (the
-- default) keeps the historical behaviour, because vdjdb-web's Scala side matches
-- `data:image/png;base64` and would find nothing to render otherwise. See `docs/dashboard.md`.
local ASSET_PREFIX = nil
local TABLE_CLASS = "ui unstackable single line celled stripped compact small table"
-- #460. knitr emits `width="672"`/`width="1152"` on these tags and pandoc's resource embedding
-- drops the attribute when it inlines the image, which is why the old pixel-width rewrite never
-- fired: it ran on the stage where the attribute no longer existed. A style cannot be dropped that
-- way, and it is responsive, which a fixed pixel width never was.
local IMG_STYLE = 'style="max-width:100%;height:auto" '

local function is_marker(block, text)
  return block.t == "Para" and pandoc.utils.stringify(block) == text
end

local function Pandoc(doc)
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

local function CodeBlock(el)
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
  if ASSET_PREFIX then
    -- `vdjdb_summary_files/figure-html/unnamed-chunk-5-1.png` -> `<prefix>unnamed-chunk-5-1.png`.
    -- knitr's path is an artifact of the render; what ships is whatever vdjdb-web serves.
    text = text:gsub('src="[^"]*/([^"/]+%.png)"',
                     function(name) return 'src="' .. ASSET_PREFIX .. name .. '"' end)
  end
  return text
end

local function RawBlock(el)
  if el.format ~= "html" then
    return el
  end
  return pandoc.RawBlock("html", patch(el.text))
end

-- knitr emits each figure as a bare `<img ...>` line, which pandoc parses as a RawInline inside a
-- Para -- NOT a RawBlock. Filtering only blocks left every image unstyled, and the structural
-- check caught it: 8 images, 0 with a style attribute.
local function RawInline(el)
  if el.format ~= "html" then
    return el
  end
  return pandoc.RawInline("html", patch(el.text))
end

-- Two passes, in this order and not by luck: pandoc walks element filters before the document-level
-- `Pandoc` function, so reading `asset-prefix` there left every `src` unrewritten while the style
-- was applied -- measured, it did. The first pass exists only to read the option.
return {
  {
    Meta = function(meta)
      if meta["asset-prefix"] then
        ASSET_PREFIX = pandoc.utils.stringify(meta["asset-prefix"])
      end
      return meta
    end,
  },
  {
    CodeBlock = CodeBlock,
    RawBlock = RawBlock,
    RawInline = RawInline,
    Pandoc = Pandoc,
  },
}
