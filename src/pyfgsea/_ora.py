"""Over-representation analysis helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import pandas as pd
from scipy.stats import hypergeom

from ._params import merge_params, pop_required
from ._utils import coerce_pathways, p_adjust_bh, prepare_pathways, validate_universe


FORA_KEYS = {"pathways", "genes", "universe", "min_size", "max_size"}


def fora(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    genes: Iterable[Any] | None = None,
    universe: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run over-representation analysis using the hypergeometric test."""

    merged = merge_params(params, allowed_keys=FORA_KEYS, context="fora", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "fora")
    genes = pop_required(merged, "genes", genes, "fora")
    universe = pop_required(merged, "universe", universe, "fora")
    min_size = int(merged.pop("min_size", 1))
    universe_list = validate_universe([str(gene) for gene in universe], "universe")
    max_size = int(merged.pop("max_size", max(len(universe_list) - 1, 0)))
    filtered, sizes, genes_by_pathway = prepare_pathways(pathways, universe_list, min_size=min_size, max_size=max_size)
    if not filtered:
        return pd.DataFrame(
            {
                "pathway": pd.Series(dtype=str),
                "pval": pd.Series(dtype=float),
                "padj": pd.Series(dtype=float),
                "fold_enrichment": pd.Series(dtype=float),
                "overlap": pd.Series(dtype=int),
                "size": pd.Series(dtype=int),
                "overlap_genes": pd.Series(dtype=object),
            }
        )

    universe_set = set(universe_list)
    gene_set = {str(gene) for gene in genes if str(gene) in universe_set}
    if len(gene_set) < len(set(str(gene) for gene in genes)):
        import warnings

        warnings.warn("Not all input genes belong to the universe; those genes were removed.", RuntimeWarning)

    rows: list[dict[str, Any]] = []
    population = len(universe_list)
    draws = len(gene_set)
    for name in filtered:
        pathway_genes = set(genes_by_pathway[name])
        overlap_genes = sorted(pathway_genes & gene_set)
        overlap = len(overlap_genes)
        size = sizes[name]
        pval = float(hypergeom.sf(overlap - 1, population, size, draws))
        expected = size / population if population else 0
        observed = overlap / draws if draws else 0
        rows.append(
            {
                "pathway": name,
                "pval": pval,
                "fold_enrichment": observed / expected if expected else float("nan"),
                "overlap": overlap,
                "size": size,
                "overlap_genes": overlap_genes,
            }
        )
    result = pd.DataFrame(rows)
    result["padj"] = p_adjust_bh(result["pval"])
    return result[["pathway", "pval", "padj", "fold_enrichment", "overlap", "size", "overlap_genes"]].sort_values(
        "pval", kind="mergesort"
    ).reset_index(drop=True)


COLLAPSE_ORA_KEYS = {"fora_res", "pathways", "genes", "universe", "pval_threshold", "overlap_threshold"}


def collapse_pathways_ora(
    fora_res: pd.DataFrame | None = None,
    pathways: Mapping[str, Iterable[Any]] | None = None,
    genes: Iterable[Any] | None = None,
    universe: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Collapse redundant ORA pathways using overlap with better-ranked pathways."""

    merged = merge_params(params, allowed_keys=COLLAPSE_ORA_KEYS, context="collapse_pathways_ora", **kwargs)
    fora_res = pop_required(merged, "fora_res", fora_res, "collapse_pathways_ora")
    pathways = pop_required(merged, "pathways", pathways, "collapse_pathways_ora")
    _ = pop_required(merged, "genes", genes, "collapse_pathways_ora")
    _ = pop_required(merged, "universe", universe, "collapse_pathways_ora")
    pval_threshold = float(merged.pop("pval_threshold", 0.05))
    overlap_threshold = float(merged.pop("overlap_threshold", 0.5))
    pathway_map = {name: set(values) for name, values in coerce_pathways(pathways).items()}
    ordered = fora_res.sort_values("pval", kind="mergesort")
    ordered = ordered[ordered["pval"] <= pval_threshold]
    parent: dict[str, str | None] = {str(name): None for name in ordered["pathway"]}
    main: list[str] = []
    for name in ordered["pathway"].astype(str):
        if parent.get(name) is not None:
            continue
        main.append(name)
        genes1 = pathway_map.get(name, set())
        for other in ordered["pathway"].astype(str):
            if other == name or parent.get(other) is not None:
                continue
            genes2 = pathway_map.get(other, set())
            denom = min(len(genes1), len(genes2)) or 1
            if len(genes1 & genes2) / denom >= overlap_threshold:
                parent[other] = name
    return {"main_pathways": main, "parent_pathways": parent}
