"""Shared data coercion and statistical helpers."""

from __future__ import annotations

import math
import warnings
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd


SCORE_TYPES = {"std", "pos", "neg"}
ENGINES = {"auto", "python", "compiled"}


def validate_score_type(score_type: str) -> str:
    if score_type not in SCORE_TYPES:
        raise ValueError("score_type must be one of: std, pos, neg.")
    return score_type


def validate_engine(engine: str) -> str:
    if engine not in ENGINES:
        raise ValueError("engine must be one of: auto, python, compiled.")
    if engine == "compiled":
        raise NotImplementedError(
            "Compiled acceleration is not built in this v1 package yet; use engine='auto' or engine='python'."
        )
    return "python"


def p_adjust_bh(pvalues: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""

    pvals = np.asarray(pvalues, dtype=float)
    adjusted = np.full(pvals.shape, np.nan, dtype=float)
    finite = np.isfinite(pvals)
    if not finite.any():
        return adjusted

    finite_p = pvals[finite]
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    n = len(ranked)
    raw = ranked * n / np.arange(1, n + 1)
    monotone = np.minimum.accumulate(raw[::-1])[::-1]
    monotone = np.clip(monotone, 0, 1)
    restored = np.empty_like(monotone)
    restored[order] = monotone
    adjusted[finite] = restored
    return adjusted


def coerce_stats(
    stats: Any,
    *,
    gene_col: str | None = None,
    value_col: str | None = None,
) -> pd.Series:
    """Coerce ranked statistics into a named numeric Series."""

    if isinstance(stats, pd.Series):
        series = stats.copy()
    elif isinstance(stats, pd.DataFrame):
        if gene_col is None:
            gene_col = "gene" if "gene" in stats.columns else stats.columns[0]
        if value_col is None:
            candidates = [c for c in ("stat", "score", "value", "t") if c in stats.columns]
            value_col = candidates[0] if candidates else stats.columns[1]
        series = pd.Series(stats[value_col].to_numpy(), index=stats[gene_col].astype(str).to_numpy(), name=value_col)
    elif isinstance(stats, Mapping):
        series = pd.Series(dict(stats), dtype=float)
    else:
        series = pd.Series(stats, dtype=float)

    if series.index is None or isinstance(series.index, pd.RangeIndex):
        raise ValueError("stats must have gene names; pass a Series, mapping, or DataFrame with gene identifiers.")
    series.index = series.index.astype(str)
    series = pd.to_numeric(series, errors="coerce")
    validate_universe(series.index, "stats index")
    if not np.isfinite(series.to_numpy(dtype=float)).all():
        raise ValueError("Not all stats values are finite numbers.")
    return series.astype(float)


def prepare_ranked_stats(
    stats: Any,
    *,
    score_type: str = "std",
    gsea_param: float = 1.0,
    gene_col: str | None = None,
    value_col: str | None = None,
) -> pd.Series:
    """Validate and sort stats for GSEA calculations."""

    score_type = validate_score_type(score_type)
    series = coerce_stats(stats, gene_col=gene_col, value_col=value_col)
    ties = int(pd.Series(series[series != 0]).duplicated().sum())
    if ties:
        warnings.warn(
            f"There are ties in the preranked stats ({round(ties * 100 / len(series), 2)}% of the list).",
            RuntimeWarning,
            stacklevel=2,
        )
    if bool((series > 0).all()) and score_type == "std":
        warnings.warn(
            'All values in stats are greater than zero and score_type is "std"; consider score_type="pos".',
            RuntimeWarning,
            stacklevel=2,
        )
    sorted_series = series.sort_values(ascending=False, kind="mergesort")
    if gsea_param < 0:
        raise ValueError("gsea_param must be non-negative.")
    return sorted_series


def coerce_expression(expression: Any) -> pd.DataFrame:
    """Coerce an expression matrix into a finite DataFrame with named genes."""

    if isinstance(expression, pd.DataFrame):
        frame = expression.copy()
    else:
        frame = pd.DataFrame(expression)
    if isinstance(frame.index, pd.RangeIndex):
        raise ValueError("expression must have gene names in its index.")
    frame.index = frame.index.astype(str)
    validate_universe(frame.index, "expression index")
    frame = frame.apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("Not all expression values are finite numbers.")
    return frame.astype(float)


def validate_universe(universe: Iterable[Any], label: str) -> list[str]:
    values = [str(value) if value is not None else None for value in universe]
    if any(value is None or value == "nan" for value in values):
        raise ValueError(f"NAs in {label} are not allowed.")
    if any(value == "" for value in values):
        raise ValueError(f"Empty strings are not allowed in {label}.")
    if len(set(values)) != len(values):
        raise ValueError(f"Duplicate values in {label} are not allowed.")
    return [str(value) for value in values]


def coerce_pathways(pathways: Mapping[str, Iterable[Any]]) -> dict[str, list[str]]:
    if not isinstance(pathways, Mapping):
        raise TypeError("pathways must be a mapping from pathway names to gene collections.")
    return {str(name): [str(gene) for gene in genes if gene is not None and str(gene) != ""] for name, genes in pathways.items()}


def prepare_pathways(
    pathways: Mapping[str, Iterable[Any]],
    universe: Sequence[str],
    *,
    min_size: int,
    max_size: int,
) -> tuple[dict[str, list[int]], dict[str, int], dict[str, list[str]]]:
    """Filter pathways and convert genes to universe positions."""

    if min_size < 1:
        min_size = 1
    max_size = min(max_size, max(len(universe) - 1, 0))
    universe = validate_universe(universe, "universe")
    positions = {gene: i for i, gene in enumerate(universe)}

    filtered: dict[str, list[int]] = {}
    sizes: dict[str, int] = {}
    genes_by_pathway: dict[str, list[str]] = {}
    for name, genes in coerce_pathways(pathways).items():
        seen: set[int] = set()
        selected: list[int] = []
        selected_genes: list[str] = []
        for gene in genes:
            if gene in positions and positions[gene] not in seen:
                seen.add(positions[gene])
                selected.append(positions[gene])
                selected_genes.append(gene)
        size = len(selected)
        if min_size <= size <= max_size:
            filtered[name] = selected
            sizes[name] = size
            genes_by_pathway[name] = selected_genes
    return filtered, sizes, genes_by_pathway


def rng_from_seed(seed: int | None) -> np.random.Generator:
    return np.random.default_rng(seed)


def safe_log2err_from_count(count: int, n_perm: int) -> float:
    if count < 0 or n_perm <= 0:
        return math.nan
    return float(math.sqrt(max(0.0, 1.0 / (count + 1) - 1.0 / (n_perm + 1))) / math.log(2))
