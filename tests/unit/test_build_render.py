# tests/unit/test_build_render.py
"""The shell every page shares, and the escaping every page depends on."""

from chd_atlas.build.render import Link, Row, chip, data_table, definition_list, document


def test_a_document_is_self_contained_and_carries_no_external_reference() -> None:
    """No `<link>`, `<script src>` or `@import`: the page renders with no
    network request beyond the one that fetched it, and a strict-CSP or offline
    reader sees the same site.
    """
    page = document(title="CHD Atlas", root="", body="<p>hi</p>")

    assert page.startswith("<!doctype html>")
    assert "<link" not in page
    assert "<script src" not in page
    assert "@import" not in page
    assert "<style>" in page


def test_the_title_is_escaped() -> None:
    assert "<script>" not in document(title="<script>x</script>", root="", body="")


def test_a_nested_page_reaches_the_root_through_the_prefix() -> None:
    """A gene page lives in `genes/`, so its links back need `../`. Passed in
    rather than computed here: `render.py` knows nothing about where a caller
    put its page.
    """
    assert 'href="../index.html"' in document(title="TBX5", root="../", body="")
    assert 'href="index.html"' in document(title="CHD Atlas", root="", body="")


def test_every_cell_of_a_table_is_escaped_including_a_link() -> None:
    """`data_table` takes text and `Link`, never HTML.

    A curator-supplied or mirror-supplied string reaching a page as markup is
    the one way this build turns data into script, and `bundles.py` and
    `search.py` already treat both as data for the same reason. Taking
    pre-rendered HTML here would move that decision to every call site.
    """
    html = data_table(
        headers=["source", "report"],
        rows=[Row(cells=("<b>clingen</b>", Link(text="open <x>", href="https://e.org/?a=1&b=2")))],
    )

    assert "<b>" not in html
    assert "&lt;b&gt;clingen&lt;/b&gt;" in html
    assert "&lt;x&gt;" in html
    assert "https://e.org/?a=1&amp;b=2" in html


def test_a_quote_in_an_href_cannot_close_the_attribute_early() -> None:
    """The escape that stops a URL from becoming a second attribute.

    `html.escape` rewrites `"` only under `quote=True`, its default; passing
    `quote=False` anywhere in `render.py` would still escape every angle bracket
    — so `test_every_cell_of_a_table_is_escaped_including_a_link` above stays
    green — while letting any href close its own attribute and open an event
    handler. Asserted against the whole anchor, because a bare `"` is legal in
    the surrounding markup and only its position inside the value is the defect.
    """
    href = 'https://e.org/?a=" onmouseover="alert(1)'
    html = data_table(headers=["report"], rows=[Row(cells=(Link(text="open", href=href),))])

    assert '<a href="https://e.org/?a=&quot; onmouseover=&quot;alert(1)">open</a>' in html
    assert ' onmouseover="' not in html


def test_a_row_attribute_name_and_value_are_both_escaped() -> None:
    """`data-*` carries curated text too — lesion group names reach it verbatim.

    The name is escaped as well as the value, because both come from the caller
    and an unescaped name is the shorter route to the same escape: no `<` is
    needed, only a `"` to close the attribute and start another one the row
    never declared.
    """
    html = data_table(
        headers=["gene"],
        rows=[Row(cells=("TBX5",), attributes=(('lesion" onclick="x', '<b> "septal"'),))],
    )

    assert ' onclick="' not in html
    assert "<b>" not in html
    assert 'data-lesion&quot; onclick=&quot;x="&lt;b&gt; &quot;septal&quot;"' in html


def test_a_table_with_no_rows_still_renders_its_header() -> None:
    """A gene with no variants and a table that was never written must not read
    the same to a reader. The header names the columns, so an empty body is an
    answer; a missing table is a defect that looks like one.
    """
    html = data_table(headers=["variant", "consequence"], rows=[])

    assert "<thead>" in html
    assert "<th>variant</th><th>consequence</th>" in html
    assert "<tbody></tbody>" in html


def test_a_chip_and_a_definition_list_escape_their_values() -> None:
    assert "<i>" not in chip(text="<i>definitive</i>", kind="definitive")
    assert "<i>" not in definition_list([("<i>k</i>", "<i>v</i>")])
