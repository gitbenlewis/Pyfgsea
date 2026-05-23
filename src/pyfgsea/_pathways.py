"""Pathway collection helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from ._params import merge_params, pop_required


GMT_KEYS = {"gmt_file"}
WRITE_GMT_KEYS = {"pathways", "gmt_file"}
REACTOME_KEYS = {"genes", "organism", "source"}
MAP_IDS_KEYS = {"x", "keys", "column", "keytype", "missing"}


def gmt_pathways(
    gmt_file: str | Path | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, list[str]]:
    """Read pathways from a GMT file."""

    merged = merge_params(params, allowed_keys=GMT_KEYS, context="gmt_pathways", **kwargs)
    gmt_file = pop_required(merged, "gmt_file", gmt_file, "gmt_pathways")
    result: dict[str, list[str]] = {}
    with Path(gmt_file).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            result[parts[0]] = [gene for gene in parts[2:] if gene]
    return result


def write_gmt_pathways(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    gmt_file: str | Path | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> None:
    """Write a pathway mapping to a GMT file."""

    merged = merge_params(params, allowed_keys=WRITE_GMT_KEYS, context="write_gmt_pathways", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "write_gmt_pathways")
    gmt_file = pop_required(merged, "gmt_file", gmt_file, "write_gmt_pathways")
    output = Path(gmt_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for name, genes in pathways.items():
            handle.write("\t".join([str(name), "NA", *[str(gene) for gene in genes]]) + "\n")


def reactome_pathways(
    genes: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, list[str]]:
    """Return Reactome pathways for a set of genes using optional gseapy support."""

    merged = merge_params(params, allowed_keys=REACTOME_KEYS, context="reactome_pathways", **kwargs)
    genes = pop_required(merged, "genes", genes, "reactome_pathways")
    organism = merged.pop("organism", "Mouse")
    source = merged.pop("source", f"Reactome_{organism}")
    try:
        import gseapy as gp  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "reactome_pathways requires optional package 'gseapy'. "
            "Install it or pass a local GMT file through gmt_pathways()."
        ) from exc
    library = gp.get_library(name=source, organism=organism)
    gene_set = {str(gene) for gene in genes}
    result: dict[str, list[str]] = {}
    for name, values in library.items():
        overlap = gene_set & {str(gene) for gene in values}
        if overlap:
            result[name] = sorted(overlap)
    return result


def map_ids_list(
    x: Mapping[Any, Any] | Any = None,
    keys: Mapping[str, Iterable[Any]] | None = None,
    column: str | None = None,
    keytype: str | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, list[Any]]:
    """Map each list of identifiers using a Python mapping or callable mapper."""

    merged = merge_params(params, allowed_keys=MAP_IDS_KEYS, context="map_ids_list", **kwargs)
    x = pop_required(merged, "x", x, "map_ids_list")
    keys = pop_required(merged, "keys", keys, "map_ids_list")
    _ = merged.pop("column", column)
    _ = merged.pop("keytype", keytype)
    missing = merged.pop("missing", None)

    def map_one(value: Any) -> Any:
        if callable(x):
            return x(value)
        return x.get(value, x.get(str(value), missing))

    return {str(name): [map_one(value) for value in values] for name, values in keys.items()}
