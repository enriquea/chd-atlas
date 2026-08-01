# tests/unit/test_build_paths.py
from typing import Final

import pytest
from pydantic import TypeAdapter

from chd_atlas.build.paths import gene_bundle_path, slug
from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    ContrastId,
    Doi,
    FunctionalId,
    HgncId,
    ModelPhenotypeId,
    ModificationId,
    PhenotypeId,
    Pmcid,
    Pmid,
    SequenceOntologyId,
    TaxonId,
    UniprotAccession,
)

# At least one ASCII value per grammar in `identifiers.py` except `Doi`, which
# slug provably cannot separate and which has its own test below, plus the two
# values closest to colliding: `HGNC:11604` against a `ContrastId` spelling the
# same characters. Each is paired with its type so the injectivity claim is
# pinned against strings the grammars actually admit rather than against
# plausible-looking inventions. ASCII because injectivity holds only there —
# `test_slug_cannot_separate_non_ascii_digits_so_injectivity_is_ascii_only`.
SAMPLE: Final = [
    (HgncId, "HGNC:11604"),
    (Pmid, "PMID:8988165"),
    (Pmcid, "PMC1234567"),
    (PhenotypeId, "HP:0001631"),
    (PhenotypeId, "MONDO:0007259"),
    (ModelPhenotypeId, "MP:0010402"),
    (ModelPhenotypeId, "ZP:0000706"),
    (SequenceOntologyId, "SO:0001587"),
    (ModificationId, "MOD:00046"),
    (TaxonId, "NCBITaxon:9606"),
    (AssertionId, "CHDA:AST:0000001"),
    (FunctionalId, "CHDA:FUN:0000042"),
    (UniprotAccession, "Q99593"),
    (UniprotAccession, "Q99593-2"),
    (AccessionId, "PXD012345"),
    (AccessionId, "GSE123456"),
    (AccessionId, "E-MTAB-1234"),
    (AccessionId, "EGAS00001000123"),
    (ContrastId, "tof_vs_control"),
    (ContrastId, "hgnc_11604"),
]


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("HGNC:11604", "HGNC_11604"),
        ("PMID:8988165", "PMID_8988165"),
        ("HP:0001631", "HP_0001631"),
        # Zygosity.NOT_APPLICABLE is the value "n/a"; an unslugged slash would
        # silently become a directory separator in any path built from it.
        ("n/a", "n_a"),
        ("PXD012345", "PXD012345"),
        # `ContrastId` is built from underscore-joined words, and slug leaves an
        # underscore as itself rather than rewriting it to something else.
        ("tof_vs_control", "tof_vs_control"),
        # The dot is safe, so a DOI keeps its "10.1000" and only the slash moves.
        # `Doi` is the only grammar containing a dot at all.
        ("10.1000/ab.cd", "10.1000_ab.cd"),
        # The hyphen is safe as well: a UniProt isoform suffix and an
        # ArrayExpress accession both survive whole. `slug`'s claim that it
        # rewrites only the colon and the underscore depends on this, so a
        # mutation dropping "-" from the safe set has to fail somewhere.
        ("Q99593-2", "Q99593-2"),
        ("E-MTAB-1234", "E-MTAB-1234"),
    ],
)
def test_slug_replaces_every_character_unsafe_in_a_path(value: str, expected: str) -> None:
    assert slug(value) == expected


def test_the_injectivity_sample_is_drawn_from_the_real_grammars() -> None:
    """`slug` claims a property of `identifiers.py`, so the sample must come from it.

    Left unchecked, a sample drifts into strings no grammar admits — "HGNC_1"
    reads like an identifier and is not one — and the test below would then pin
    a property of values that cannot occur.
    """
    for type_, value in SAMPLE:
        assert TypeAdapter(type_).validate_python(value) == value


def test_slug_is_injective_over_the_ascii_values_of_every_grammar_but_doi() -> None:
    """Two distinct ASCII identifiers must not collapse to one filename.

    Within ASCII the only characters slug rewrites are the colon and the
    underscore, and no two values differ solely by swapping the two: the
    underscore appears only in `ContrastId`, which admits no capital, while every
    grammar bearing a colon starts with one — and the hyphen and dot, which the
    parametrised cases above pin, are left alone entirely.

    Two things this does not cover, each with its own test below: `Doi`, and any
    non-ASCII value of any grammar.
    """
    values = [value for _, value in SAMPLE]
    assert len({slug(value) for value in values}) == len(values)


def test_slug_cannot_separate_dois_so_no_path_may_be_built_from_one() -> None:
    r"""The one grammar excluded outright, pinned so the exclusion stays true.

    `Doi` is `^10\.\d{4,9}/\S+$`: the suffix admits every non-whitespace
    character, so slash, colon and underscore all arrive here and all leave as
    "_". This fails if something later relies on slug to name a DOI uniquely.
    """
    colliding = ["10.1000/ab/cd", "10.1000/ab:cd", "10.1000/ab_cd"]
    for value in colliding:
        assert TypeAdapter(Doi).validate_python(value) == value
    assert len(set(colliding)) == 3
    assert len({slug(value) for value in colliding}) == 1


def test_slug_cannot_separate_non_ascii_digits_so_injectivity_is_ascii_only() -> None:
    r"""The second exception, and the one that reaches `HgncId` itself.

    `\d` is Unicode-aware, so `^HGNC:\d+$` accepts far more than "HGNC:" and
    ASCII digits. The two values below are both accepted by the project's own
    validator, are genuinely distinct, and share one bundle path — which makes
    the injectivity above an ASCII property rather than a total one.

    This costs no gene page today: no curated source emits a non-ASCII digit,
    and `Emitter._write` raises on the second write to a path rather than
    replacing the first, so a collision would fail the build. That is a backstop
    elsewhere, not a property of `slug`, which is why it is pinned here.
    """
    # ARABIC-INDIC DIGIT ONE and TWO (U+0661, U+0662). They are near-identical to
    # ASCII in some fonts; they are not typos and must not be "corrected".
    colliding = ["HGNC:1160١", "HGNC:1160٢"]
    for value in colliding:
        assert TypeAdapter(HgncId).validate_python(value) == value
    assert len(set(colliding)) == 2
    assert len({slug(value) for value in colliding}) == 1
    assert len({gene_bundle_path(value) for value in colliding}) == 1


def test_gene_bundle_path_is_relative_and_forward_slashed() -> None:
    assert gene_bundle_path("HGNC:11604") == "genes/HGNC_11604.json"
