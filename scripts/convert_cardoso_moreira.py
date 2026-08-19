#!/usr/bin/env python3
"""Build the bulk developmental expression mirrors from the Expression Atlas.

Source: **E-MTAB-6814**, the EMBL-EBI Expression Atlas reprocessing of
Cardoso-Moreira et al. 2019 (*Nature*) -- human developmental transcriptomes
across seven organs.

**The published values are the Atlas's, not the paper's, and that distinction
is load-bearing.** `E-MTAB-6814-analysis-methods.tsv` records the Atlas
pipeline: iRAP 1.0.1, HISAT2 2.1.0 against Ensembl release 107, FeatureCounts
1.6.2, then limma quantile normalisation within each set of biological
replicates. The paper's own supplementary tables report RPKM from a different
pipeline. The two are not interchangeable and the provenance must name which
produced a number. (A stale note survives in this repository's history: the
commit adding `rpkm` to the unit vocabulary says the source "publishes RPKM" --
true of the paper, false of the Atlas mirror this script reads. `tpm` is what
is emitted here.)

Why TPM rather than FPKM, both of which the Atlas publishes: TPM sums to a
constant per sample, so it is the more defensible quantity to compare across
samples. FPKM's across-gene sum varies between samples.

## What the Atlas gives us for free

Each matrix cell is not a scalar but a five-value quantile summary of that
gene's assay group -- ``min,q25,median,q75,max`` (e.g. ``8,15,18,20,21``). So
`median_abundance`, `q25` and `q75` are read off directly rather than computed,
and `configuration.xml` maps each assay-group id to a
``"<stage>; <organ>"`` label plus the assays it aggregates, which yields
`tissue`, `stage` and `n_samples`.

## The two outputs, and why they are different populations

- **profiles** -- only the genes in `mirrors/genes.tsv` (the registry, 154
  today). Deliberately wider than the publication gate: widening D21 later must
  not require re-mirroring.
- **profile_quantiles** -- breakpoints over the matrix's **entire** gene
  universe (58,735), because a percentile's denominator is the whole
  transcriptome. Restricting it to the registry would publish a rank against a
  population no reader means.

Gene identity is joined on `ensembl_gene`, never on a symbol. Measured
2026-08-15: all 154 registry genes are present in the matrix, none missing.
Any that ever go missing are named on stdout rather than dropped silently.
"""

from __future__ import annotations

import argparse
import csv
import sys
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Sequence
from pathlib import Path

# The quantile estimator, named because D39(a) requires the method to travel
# with the number, and implemented here rather than delegated to a library so
# the definition is auditable from this file alone. `linear` is the R type-7 /
# numpy-default rule: with the values sorted ascending and h = (n - 1) * p,
# the estimate is a linear interpolation between the values at floor(h) and
# ceil(h). Pinned against worked literals in the unit tests -- a reproducibility
# script must not inherit a quantile convention from a dependency's version.
QUANTILE_ESTIMATOR = "linear"
PERCENTILES = list(range(101))
UNIT = "tpm"


class AssayGroup:
    """One assay group: a (stage, tissue) cell and the assays it aggregates."""

    __slots__ = ("group_id", "stage", "tissue", "n_samples")

    def __init__(self, group_id: str, stage: str, tissue: str, n_samples: int) -> None:
        self.group_id = group_id
        self.stage = stage
        self.tissue = tissue
        self.n_samples = n_samples


def parse_assay_groups(config_xml: Path) -> dict[str, AssayGroup]:
    """Map each assay-group id to its (stage, tissue) and replicate count.

    Labels read ``"10 week post conception; heart"``. Split on the first
    semicolon only: an organ name never contains one, and a stage token might.
    A label without a semicolon is skipped and named -- it is a group the
    experiment defines on some axis other than stage-and-organ, and guessing
    which half is which would silently mislabel a whole column.
    """
    groups: dict[str, AssayGroup] = {}
    root = ET.parse(config_xml).getroot()
    for element in root.iter("assay_group"):
        group_id = element.get("id")
        label = element.get("label")
        if not group_id or not label:
            continue
        if ";" not in label:
            print(f"skipped group {group_id}: label has no '; ' split: {label!r}", file=sys.stderr)
            continue
        stage, tissue = (part.strip() for part in label.split(";", 1))
        n_samples = sum(1 for _ in element.iter("assay"))
        groups[group_id] = AssayGroup(group_id, stage, tissue, n_samples)
    return groups


def parse_cell(cell: str) -> tuple[float, float, float] | None:
    """(q25, median, q75) from the Atlas's ``min,q25,median,q75,max`` quintuple.

    Returns None for an absent measurement rather than zero: a gene the Atlas
    did not quantify in a group is not a gene measured at zero, and collapsing
    the two would publish a floor-level abundance the source never reported.

    A cell with a single value is accepted as a degenerate quintuple (all five
    equal), which is what the Atlas emits for a one-assay group.
    """
    cell = cell.strip()
    if not cell or cell.upper() in {"NA", "NAN"}:
        return None
    parts = cell.split(",")
    try:
        values = [float(part) for part in parts]
    except ValueError:
        return None
    if len(values) == 1:
        return (values[0], values[0], values[0])
    if len(values) != 5:
        return None
    return (values[1], values[2], values[3])


def profile_row(
    dataset: str, gene: str, group: AssayGroup, q25: float, median: float, q75: float
) -> dict[str, object]:
    """One `mirrors/profiles/` row.

    `q25`/`q75` are emitted as null below three replicates: quartiles of two
    points are not quartiles, and a median of two is a mean of two with none of
    the robustness the word implies.
    """
    has_quartiles = group.n_samples >= 3
    return {
        "dataset": dataset,
        "gene": gene,
        "tissue": group.tissue,
        "stage": group.stage,
        "median_abundance": median,
        "unit": UNIT,
        "q25": q25 if has_quartiles else None,
        "q75": q75 if has_quartiles else None,
        "n_samples": group.n_samples,
    }


def quantile_linear(sorted_values: Sequence[float], percentile: int) -> float:
    """The R type-7 / numpy-default quantile of an already-sorted sequence.

    Written out rather than imported so the estimator this atlas publishes is
    readable in the file that produced the numbers. `sorted_values` must be
    ascending; sorting inside would hide an O(n log n) cost per percentile.
    """
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    position = (n - 1) * (percentile / 100.0)
    lower = int(position)
    upper = min(lower + 1, n - 1)
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def quantile_rows(
    dataset: str, group: AssayGroup, medians: list[float]
) -> Iterator[dict[str, object]]:
    """101 breakpoints for one (tissue, stage) cell, over the whole universe.

    Computed on every gene the matrix quantified in this group, not on the
    registry subset -- the denominator of a published percentile is the source's
    own gene universe, and `n_genes` states it.
    """
    ordered = sorted(medians)
    n_genes = len(ordered)
    for percentile in PERCENTILES:
        value = quantile_linear(ordered, percentile)
        yield {
            "dataset": dataset,
            "tissue": group.tissue,
            "stage": group.stage,
            "percentile": percentile,
            "value": round(float(value), 6),
            "unit": UNIT,
            "n_genes": n_genes,
        }


def _write(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    """Write a TSV with the empty string for null, matching `read_table`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
            restval="",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out-profiles", type=Path, required=True)
    parser.add_argument("--out-quantiles", type=Path, required=True)
    args = parser.parse_args()

    groups = parse_assay_groups(args.config)
    print(f"assay groups: {len(groups)}", file=sys.stderr)

    ensembl_to_hgnc: dict[str, str] = {}
    with args.registry.open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            if record.get("ensembl_gene"):
                ensembl_to_hgnc[record["ensembl_gene"]] = record["hgnc_id"]
    print(f"registry genes with an ensembl id: {len(ensembl_to_hgnc)}", file=sys.stderr)

    profiles: list[dict[str, object]] = []
    collected: dict[str, list[float]] = {gid: [] for gid in groups}
    seen: set[str] = set()

    with args.matrix.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        columns = {index: header[index] for index in range(2, len(header))}
        for record in reader:
            ensembl = record[0]
            hgnc = ensembl_to_hgnc.get(ensembl)
            if hgnc is not None:
                seen.add(ensembl)
            for index, group_id in columns.items():
                group = groups.get(group_id)
                if group is None or index >= len(record):
                    continue
                parsed = parse_cell(record[index])
                if parsed is None:
                    continue
                q25, median, q75 = parsed
                collected[group_id].append(median)
                if hgnc is not None:
                    profiles.append(profile_row(args.dataset, hgnc, group, q25, median, q75))

    missing = sorted(set(ensembl_to_hgnc) - seen)
    if missing:
        print(f"WARNING: {len(missing)} registry genes absent from the matrix:", file=sys.stderr)
        for ensembl in missing:
            print(f"  {ensembl} ({ensembl_to_hgnc[ensembl]})", file=sys.stderr)
    else:
        print("every registry gene is present in the matrix", file=sys.stderr)

    quantiles: list[dict[str, object]] = []
    for group_id, group in groups.items():
        medians = collected[group_id]
        if not medians:
            print(f"WARNING: group {group_id} quantified no gene; no grid", file=sys.stderr)
            continue
        quantiles.extend(quantile_rows(args.dataset, group, medians))

    # Sorted to each table's declared sort key: SORT001 is an error, not a warning.
    profiles.sort(key=lambda r: (r["dataset"], r["gene"], r["tissue"], r["stage"]))
    quantiles.sort(key=lambda r: (r["dataset"], r["tissue"], r["stage"], r["percentile"]))

    _write(
        args.out_profiles,
        [
            "dataset",
            "gene",
            "tissue",
            "stage",
            "median_abundance",
            "unit",
            "q25",
            "q75",
            "n_samples",
        ],
        profiles,
    )
    _write(
        args.out_quantiles,
        ["dataset", "tissue", "stage", "percentile", "value", "unit", "n_genes"],
        quantiles,
    )

    print(f"profiles rows:  {len(profiles)}", file=sys.stderr)
    print(f"quantile rows:  {len(quantiles)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
