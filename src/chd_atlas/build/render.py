# src/chd_atlas/build/render.py
"""The document shell and the HTML primitives every page is built from.

This module knows nothing about genes, assertions or validity. It knows how to
turn text into a page. It holds every *primitive* that escapes -- so escaping is
decided here rather than at each call site -- but it is **not** the only place
that produces markup: `pages.py` and `landing.py` both assemble their own
elements around these primitives, and this paragraph claimed otherwise until
review on #15 measured it. The claim mattered, because it is what made an
unescaped `<option value>` in `pages.py` look impossible. What is true is
narrower and is the rule to hold: **a value derived from curated or mirrored
data must reach markup only through a primitive here.**

Self-contained by construction, like every artifact this build writes: no
`<link>`, no `<script src>`, no `@import`, no `url(` in the stylesheet, no
external font, no `<img>` -- only an inline `<style>` and, on the browse page,
an inline `<script>`, both literals in this file. A page renders with no network
request beyond the one that fetched it, from a local checkout and from the
`dist/` directory a reviewer unpacks.

That is not the same as surviving a strict Content-Security-Policy, which this
paragraph also used to claim. Measured on a real build: `default-src 'self'`
blocks an inline `<style>` and an inline `<script>` exactly as it blocks a
remote one, so under such a policy the browse filter and all styling die. GitHub
Pages sets no CSP, so nothing breaks on deploy; a consumer self-hosting these
pages behind one needs a hash or nonce allowance.

Deterministic for the reason every other builder here is: nothing reads a clock,
and every value interpolated comes from the caller.

`data_table`, `chip` and `definition_list` take **text**, never HTML. A curated
gene symbol or a mirrored disease label reaching a page as markup is the one way
this build turns data into script -- `bundles.py` and `search.py` already treat
curator- and mirror-supplied strings as data for the same reason -- and taking
pre-rendered HTML would move that decision to every caller. `Link` exists so a
cell can carry a URL without any caller needing to write an `<a>` tag.

`document` renders `RESEARCH_USE_NOTICE` itself, so the statement is on the page
before any caller's body is. It used to live in `landing.py`, which meant it
reached `index.html` and nothing else: measured 2026-08-04 against a real build,
`grep -l "not a clinical decision-support tool" ` matched 1 of the 25 HTML files
the build writes, and the 24 that carry a gene-level classification -- 23 gene
pages and the browse page -- carried no research-use statement at all. A
clinician arriving at `genes/HGNC_4173.html` from a search engine met `GATA4`, an
HGNC id and a green `definitive` chip, first in DOM order and filling the first
screen under 46rem. Putting it here is what makes it structurally impossible to
add a page kind without it.
"""

from __future__ import annotations

import html
import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

# What `data_table` accepts as a `data-*` attribute name. Deliberately narrower
# than HTML permits: every name this build passes is a literal, so the only way
# to reach a rejection is for a caller to have derived one from data.
_ATTRIBUTE_NAME: Final = re.compile(r"[a-z][a-z0-9-]*")

# Moved here from `landing.py`, which now renders through `document` below.
# One stylesheet rather than one per page kind: three `<style>` blocks are three
# things that drift, and the landing page's palette is the one this extends.
STYLESHEET: Final = """
  :root {
    color-scheme: light dark;
    --fg: #1a1a1a; --bg: #ffffff; --muted: #555555; --border: #d0d0d0;
    --notice-bg: #fff4e5; --notice-border: #b45309; --link: #0b5fa5;
    --chip-bg: rgba(127, 127, 127, 0.15); --definitive: #1a7f37; --warn: #b02a37;
  }
  @media (prefers-color-scheme: dark) {
    :root {
      --fg: #e8e8e8; --bg: #14161a; --muted: #a8a8a8; --border: #3a3d42;
      --notice-bg: #3a2a12; --notice-border: #d99a3d; --link: #6cb6f5;
      --definitive: #3fb950; --warn: #f85149;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; padding: 2rem 1.25rem 4rem; background: var(--bg); color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica,
      Arial, sans-serif;
    line-height: 1.5;
  }
  main { max-width: 58rem; margin: 0 auto; }
  nav { max-width: 58rem; margin: 0 auto 1.5rem; font-size: 0.9rem; color: var(--muted); }
  h1 { margin-bottom: 0.25rem; }
  .tagline { color: var(--muted); margin-top: 0; }
  h2 { margin-top: 2.25rem; border-bottom: 1px solid var(--border); padding-bottom: 0.25rem; }
  .notice {
    max-width: 58rem; margin: 0 auto 1.5rem;
    border: 1px solid var(--notice-border); background: var(--notice-bg);
    border-radius: 6px; padding: 0.65rem 0.9rem; font-size: 0.9rem;
  }
  .notice p { margin: 0; }
  .scope-rule { color: var(--muted); font-size: 0.9rem; }
  dl { display: grid; grid-template-columns: auto 1fr; gap: 0.35rem 1rem; margin: 0; }
  dt { color: var(--muted); }
  dd { margin: 0; overflow-wrap: break-word; }
  a { color: var(--link); }
  code { background: var(--chip-bg); padding: 0.1rem 0.3rem; border-radius: 3px; }
  ul.links { list-style: none; padding: 0; margin: 0; }
  ul.links li { margin: 0.35rem 0; }
  footer { margin-top: 3rem; color: var(--muted); font-size: 0.9rem; }
  table { border-collapse: collapse; width: 100%; font-size: 0.9rem; }
  th, td { text-align: left; padding: 0.35rem 0.6rem; border-bottom: 1px solid var(--border); }
  th { color: var(--muted); font-weight: 600; }
  .chip {
    display: inline-block; background: var(--chip-bg); border-radius: 3px;
    padding: 0.1rem 0.45rem; font-size: 0.85rem; margin: 0 0.25rem 0.25rem 0;
  }
  .chip-definitive { background: var(--definitive); color: #fff; }
  .chip-warn { background: var(--warn); color: #fff; }
  .layout { display: grid; grid-template-columns: 15rem 1fr; gap: 2rem; }
  @media (max-width: 46rem) { .layout { grid-template-columns: 1fr; } }
  .rail { font-size: 0.9rem; }
  .rail h1 { font-size: 1.6rem; margin-top: 0; }
  .filters { display: flex; flex-wrap: wrap; gap: 0.5rem; margin: 1rem 0; }
  .filters input, .filters select {
    padding: 0.3rem 0.5rem; border: 1px solid var(--border); border-radius: 4px;
    background: var(--bg); color: var(--fg);
  }
"""

# The one piece of client-side behaviour on the site, and a literal so it does
# not affect determinism. It only *hides* rows that the build already rendered,
# so the browse page is complete with JavaScript disabled -- which is why the
# table is server-rendered rather than fetched (D29).
FILTER_SCRIPT: Final = """
(function () {
  var form = document.getElementById('filters');
  if (!form) return;
  var rows = Array.prototype.slice.call(
    document.querySelectorAll('#gene-table tbody tr')
  );
  var count = document.getElementById('shown');
  function apply() {
    var text = (form.querySelector('[name=q]').value || '').toLowerCase();
    var facets = Array.prototype.slice.call(form.querySelectorAll('select'));
    var shown = 0;
    rows.forEach(function (row) {
      var ok = row.getAttribute('data-search').indexOf(text) !== -1;
      facets.forEach(function (select) {
        var want = select.value;
        if (ok && want) {
          ok = (row.getAttribute('data-' + select.name) || '').split(' ').indexOf(want) !== -1;
        }
      });
      row.hidden = !ok;
      if (ok) shown += 1;
    });
    if (count) count.textContent = String(shown);
  }
  form.addEventListener('input', apply);
  apply();
})();
"""


# The project owner's own words, and the one string in this build that must
# reach a reader character for character. Kept as a single line with single
# spaces so a check for it is a plain substring check on the page's bytes rather
# than one that has to normalise a source wrap away first; `<strong>` is markup
# from this module, in the same class as `FILTER_SCRIPT` above, never a value
# from a caller.
#
# Do not reword, shorten or split this sentence. It is the site's only statement
# about what the atlas may and may not be used for, and every page now carries
# it, so a change here changes 25 pages at once.
RESEARCH_USE_NOTICE: Final = (
    "This atlas is under active development. It is "
    "<strong>not a clinical decision-support tool</strong> and must not be used to make "
    "or guide a diagnostic, treatment or any other clinical decision."
)


@dataclass(frozen=True)
class Link:
    """A table cell that is a hyperlink. `text` and `href` are both escaped."""

    text: str
    href: str


Cell = str | Link


@dataclass(frozen=True)
class Row:
    """One table row: its cells, and any `data-*` attributes it carries.

    `attributes` exists for the browse page, whose inline filter reads
    `data-lesion`, `data-confidence` and the rest off each `<tr>` rather than
    re-deriving them from the cell text. Names and values are escaped here, so
    no caller writes an attribute itself.
    """

    cells: tuple[Cell, ...]
    attributes: tuple[tuple[str, str], ...] = ()


def _cell(value: Cell) -> str:
    if isinstance(value, Link):
        return f'<a href="{html.escape(value.href)}">{html.escape(value.text)}</a>'
    return html.escape(value)


def chip(text: str, kind: str = "") -> str:
    """A small labelled pill. `kind` selects a `.chip-*` modifier class."""
    suffix = f" chip-{html.escape(kind)}" if kind else ""
    return f'<span class="chip{suffix}">{html.escape(text)}</span>'


def definition_list(items: Sequence[tuple[str, str]]) -> str:
    """A `<dl>` of label/value pairs, both escaped."""
    body = "".join(
        f"<dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd>" for label, value in items
    )
    return f"<dl>{body}</dl>"


def data_table(headers: Sequence[str], rows: Sequence[Row], table_id: str = "") -> str:
    """A `<table>` whose every cell and attribute is escaped by this function.

    An empty `rows` still renders the header, because a table with no rows and
    a table that was never written read the same to a caller otherwise.

    An attribute *name* is refused rather than escaped, because escaping cannot
    make it safe. `html.escape` rewrites `<`, `>`, `&`, `"` and `'` and leaves
    the space and the `=` alone, and a name sits outside any quote delimiter --
    so `("x onmouseover=alert(1) y", "v")` renders
    `<tr data-x onmouseover=alert(1) y="v">`, which the HTML tokenizer reads as
    three attributes, the second a live event handler, from a payload carrying
    no quote and no angle bracket. Measured with `html.parser` on 2026-08-04:
    `[('data-x', None), ('onmouseover', 'alert(1)'), ('y', 'v')]`.

    Every name the build passes is a literal in `pages.py`, so this is a guard
    on a bypassed gate in the same idiom as `encode_json`'s `allow_nan=False`:
    reaching it means a caller derived an attribute name from data, and that
    must fail rather than publish. `raise`, never `assert` -- `-O` strips
    `assert`.
    """
    ident = f' id="{html.escape(table_id)}"' if table_id else ""
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body = ""
    for row in rows:
        for name, _ in row.attributes:
            if not _ATTRIBUTE_NAME.fullmatch(name):
                raise ValueError(
                    f"data attribute name {name!r} is not [a-z][a-z0-9-]*; escaping cannot "
                    f"make a name safe, because a space in one opens a second attribute"
                )
        attributes = "".join(
            f' data-{name}="{html.escape(value)}"' for name, value in row.attributes
        )
        cells = "".join(f"<td>{_cell(cell)}</td>" for cell in row.cells)
        body += f"<tr{attributes}>{cells}</tr>"
    return f"<table{ident}><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def document(*, title: str, root: str, body: str, script: str = "") -> str:
    """The shell every page shares.

    `root` is the relative prefix from this page back to the site root -- "" for
    a page at the top level and "../" for one in `genes/`. Passed in rather than
    derived from a path, because this module has no opinion about where a caller
    put its output.

    `body` is HTML the caller has already assembled from the primitives above;
    it is the one argument not escaped here, which is why every primitive exists.
    `script` is likewise inlined verbatim and must be a literal from this module.

    `RESEARCH_USE_NOTICE` is emitted **before** `<main>`, so it precedes whatever
    the caller assembled in DOM order, in reading order and in the order a screen
    reader reaches it -- on a gene page that means it comes before the symbol,
    the HGNC id and the `definitive` chip rather than after them. It carries
    `role="note"` and no heading of its own: a gene page's `<h1>` lives inside
    `pages._rail`, and a heading here would put an `<h2>` above every `<h1>` on
    the site.
    """
    tail = f"<script>{script}</script>" if script else ""
    # Assembled ahead of the template only to keep the line under 100 characters;
    # the markup is one line in the output, exactly as it reads here.
    home = html.escape(root)
    nav = (
        f'<nav><a href="{home}index.html">CHD Atlas</a> &middot; '
        f'<a href="{home}genes/index.html">Genes</a></nav>'
    )
    notice = f'<div class="notice" role="note"><p>{RESEARCH_USE_NOTICE}</p></div>'
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{STYLESHEET}</style>
</head>
<body>
{nav}
{notice}
<main>
{body}
</main>
{tail}
</body>
</html>
"""
