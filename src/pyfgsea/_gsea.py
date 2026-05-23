"""Preranked gene set enrichment analysis."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any
import math
import warnings

import numpy as np
import pandas as pd
from scipy import special

from ._params import merge_params, pop_required
from ._utils import (
    p_adjust_bh,
    prepare_pathways,
    prepare_ranked_stats,
    rng_from_seed,
    safe_log2err_from_count,
    validate_engine,
    validate_score_type,
)


def _calc_gsea_stat_positions(
    values: np.ndarray,
    selected: Iterable[int],
    *,
    gsea_param: float = 1.0,
    return_all_extremes: bool = False,
    return_leading_edge: bool = False,
    score_type: str = "std",
) -> float | dict[str, Any]:
    score_type = validate_score_type(score_type)
    selected_positions = np.array(sorted(set(int(x) for x in selected)), dtype=int)
    if selected_positions.size == 0:
        raise ValueError("selected_stats must contain at least one gene present in stats.")
    if selected_positions[0] < 0 or selected_positions[-1] >= len(values):
        raise IndexError("selected_stats contains a position outside stats.")
    m = int(selected_positions.size)
    n = int(len(values))
    if m == n:
        raise ValueError("GSEA statistic is not defined when all genes are selected.")

    r_adj = np.abs(values[selected_positions]) ** gsea_param
    norm = float(r_adj.sum())
    if norm == 0:
        running_hit = np.arange(1, m + 1, dtype=float) / m
    else:
        running_hit = np.cumsum(r_adj) / norm

    misses_before_hit = selected_positions - np.arange(m)
    tops = running_hit - misses_before_hit / (n - m)
    if norm == 0:
        bottoms = tops - 1 / m
    else:
        bottoms = tops - r_adj / norm

    max_p = float(np.max(tops))
    min_p = float(np.min(bottoms))
    if score_type == "std":
        score = 0.0 if max_p == -min_p else (max_p if max_p > -min_p else min_p)
    elif score_type == "pos":
        score = max_p
    else:
        score = min_p

    if not return_all_extremes and not return_leading_edge:
        return float(score)

    result: dict[str, Any] = {"res": float(score)}
    if return_all_extremes:
        result["tops"] = tops
        result["bottoms"] = bottoms
    if return_leading_edge:
        if score_type == "std":
            if max_p > -min_p:
                edge = selected_positions[: np.argmax(tops) + 1]
            elif max_p < -min_p:
                edge = selected_positions[np.argmin(bottoms) :][::-1]
            else:
                edge = np.array([], dtype=int)
        elif score_type == "pos":
            edge = selected_positions[: np.argmax(tops) + 1]
        else:
            edge = selected_positions[np.argmin(bottoms) :][::-1]
        result["leading_edge_positions"] = edge.tolist()
    return result


CALC_GSEA_STAT_KEYS = {
    "stats",
    "selected_stats",
    "gsea_param",
    "return_all_extremes",
    "return_leading_edge",
    "score_type",
    "stats_gene_col",
    "stats_value_col",
}


def calc_gsea_stat(
    stats: Any = None,
    selected_stats: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> float | dict[str, Any]:
    """Calculate the enrichment score for one gene set."""

    merged = merge_params(params, allowed_keys=CALC_GSEA_STAT_KEYS, context="calc_gsea_stat", **kwargs)
    stats = pop_required(merged, "stats", stats, "calc_gsea_stat")
    selected_stats = pop_required(merged, "selected_stats", selected_stats, "calc_gsea_stat")
    gsea_param = float(merged.pop("gsea_param", 1.0))
    return_all_extremes = bool(merged.pop("return_all_extremes", False))
    return_leading_edge = bool(merged.pop("return_leading_edge", False))
    score_type = merged.pop("score_type", "std")
    stats_gene_col = merged.pop("stats_gene_col", None)
    stats_value_col = merged.pop("stats_value_col", None)

    series = prepare_ranked_stats(
        stats,
        score_type=score_type,
        gsea_param=gsea_param,
        gene_col=stats_gene_col,
        value_col=stats_value_col,
    )
    index = list(series.index)
    positions: list[int] = []
    for item in selected_stats:
        if isinstance(item, (int, np.integer)):
            positions.append(int(item))
        elif str(item) in series.index:
            positions.append(index.index(str(item)))
    result = _calc_gsea_stat_positions(
        series.to_numpy(),
        positions,
        gsea_param=gsea_param,
        return_all_extremes=return_all_extremes,
        return_leading_edge=return_leading_edge,
        score_type=score_type,
    )
    if isinstance(result, dict) and "leading_edge_positions" in result:
        result["leading_edge"] = [index[i] for i in result.pop("leading_edge_positions")]
    return result


GSEA_COMMON_KEYS = {
    "pathways",
    "stats",
    "min_size",
    "max_size",
    "score_type",
    "gsea_param",
    "seed",
    "stats_gene_col",
    "stats_value_col",
    "engine",
}

FGSEA_SIMPLE_KEYS = GSEA_COMMON_KEYS | {"n_perm", "n_jobs"}
FGSEA_MULTILEVEL_KEYS = GSEA_COMMON_KEYS | {"sample_size", "eps", "n_jobs", "n_perm_simple"}
FGSEA_KEYS = FGSEA_MULTILEVEL_KEYS | {"n_perm"}


def _empty_gsea_result(simple: bool) -> pd.DataFrame:
    columns = ["pathway", "pval", "padj"]
    if not simple:
        columns.append("log2err")
    columns.extend(["es", "nes"])
    if simple:
        columns.append("n_more_extreme")
    columns.extend(["size", "leading_edge"])
    return pd.DataFrame({column: pd.Series(dtype=object) for column in columns})


def _pathway_scores(
    pathways: Mapping[str, list[int]],
    series: pd.Series,
    *,
    gsea_param: float,
    score_type: str,
) -> tuple[dict[str, float], dict[str, list[str]]]:
    values = series.to_numpy()
    genes = list(series.index)
    scores: dict[str, float] = {}
    leading_edges: dict[str, list[str]] = {}
    for name, selected in pathways.items():
        stat = _calc_gsea_stat_positions(
            values,
            selected,
            gsea_param=gsea_param,
            return_leading_edge=True,
            score_type=score_type,
        )
        assert isinstance(stat, dict)
        scores[name] = float(stat["res"])
        leading_edges[name] = [genes[i] for i in stat["leading_edge_positions"]]
    return scores, leading_edges


def fgsea_simple(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run simple permutation-based preranked GSEA."""

    merged = merge_params(params, allowed_keys=FGSEA_SIMPLE_KEYS, context="fgsea_simple", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "fgsea_simple")
    stats = pop_required(merged, "stats", stats, "fgsea_simple")
    n_perm = int(merged.pop("n_perm", 1000))
    min_size = int(merged.pop("min_size", 1))
    max_size = int(merged.pop("max_size", 10**12))
    score_type = validate_score_type(merged.pop("score_type", "std"))
    n_jobs = int(merged.pop("n_jobs", 0))
    gsea_param = float(merged.pop("gsea_param", 1.0))
    seed = merged.pop("seed", None)
    stats_gene_col = merged.pop("stats_gene_col", None)
    stats_value_col = merged.pop("stats_value_col", None)
    engine = validate_engine(merged.pop("engine", "auto"))
    if engine != "python":
        raise AssertionError("validate_engine should normalize to python.")
    if n_jobs not in (0, 1):
        warnings.warn("n_jobs is accepted for API compatibility; the Python engine currently runs serially.", RuntimeWarning)
    if n_perm < 1:
        raise ValueError("n_perm must be at least 1.")

    series = prepare_ranked_stats(
        stats,
        score_type=score_type,
        gsea_param=gsea_param,
        gene_col=stats_gene_col,
        value_col=stats_value_col,
    )
    filtered, sizes, _ = prepare_pathways(pathways, list(series.index), min_size=min_size, max_size=max_size)
    if not filtered:
        return _empty_gsea_result(simple=True)

    actual_scores, leading_edges = _pathway_scores(filtered, series, gsea_param=gsea_param, score_type=score_type)
    values = series.to_numpy()
    rng = rng_from_seed(seed)
    rows: list[dict[str, Any]] = []
    for name, selected in filtered.items():
        size = sizes[name]
        actual = actual_scores[name]
        random_scores = np.empty(n_perm, dtype=float)
        for i in range(n_perm):
            random_selected = rng.choice(len(values), size=size, replace=False)
            random_scores[i] = float(
                _calc_gsea_stat_positions(values, random_selected, gsea_param=gsea_param, score_type=score_type)
            )
        if score_type == "pos":
            extreme = random_scores >= actual
            mode = random_scores >= 0
            denom_values = random_scores[mode]
            nes_denom = float(np.mean(denom_values)) if denom_values.size else math.nan
            nes = actual / nes_denom if nes_denom else math.nan
        elif score_type == "neg":
            extreme = random_scores <= actual
            mode = random_scores <= 0
            denom_values = random_scores[mode]
            nes_denom = abs(float(np.mean(denom_values))) if denom_values.size else math.nan
            nes = actual / nes_denom if nes_denom else math.nan
        elif actual >= 0:
            extreme = random_scores >= actual
            mode = random_scores >= 0
            denom_values = random_scores[mode]
            nes_denom = float(np.mean(denom_values)) if denom_values.size else math.nan
            nes = actual / nes_denom if nes_denom else math.nan
        else:
            extreme = random_scores <= actual
            mode = random_scores <= 0
            denom_values = random_scores[mode]
            nes_denom = abs(float(np.mean(denom_values))) if denom_values.size else math.nan
            nes = actual / nes_denom if nes_denom else math.nan
        n_more_extreme = int(np.sum(extreme))
        pval = (n_more_extreme + 1) / (n_perm + 1)
        rows.append(
            {
                "pathway": name,
                "pval": pval,
                "es": actual,
                "nes": nes,
                "n_more_extreme": n_more_extreme,
                "size": size,
                "leading_edge": leading_edges[name],
            }
        )

    result = pd.DataFrame(rows)
    result["padj"] = p_adjust_bh(result["pval"])
    result = result[["pathway", "pval", "padj", "es", "nes", "n_more_extreme", "size", "leading_edge"]]
    return result.sort_values("pval", kind="mergesort").reset_index(drop=True)


def multilevel_error(
    pval: float | np.ndarray | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> float | np.ndarray:
    """Expected log2 p-value error used by the multilevel estimator."""

    merged = merge_params(params, allowed_keys={"pval", "sample_size"}, context="multilevel_error", **kwargs)
    pval = pop_required(merged, "pval", pval, "multilevel_error")
    sample_size = int(merged.pop("sample_size", 101))
    p = np.asarray(pval, dtype=float)
    p = np.clip(p, np.finfo(float).tiny, 1)
    level = np.floor(-np.log2(p) + 1)
    single = special.polygamma(1, (sample_size + 1) / 2) - special.polygamma(1, sample_size + 1)
    err = np.sqrt(level * single) / np.log(2)
    return float(err) if np.isscalar(pval) else err


def fgsea_multilevel(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run the v1 multilevel-compatible GSEA API using the Python estimator."""

    merged = merge_params(params, allowed_keys=FGSEA_MULTILEVEL_KEYS, context="fgsea_multilevel", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "fgsea_multilevel")
    stats = pop_required(merged, "stats", stats, "fgsea_multilevel")
    sample_size = int(merged.pop("sample_size", 101))
    eps = float(merged.pop("eps", 1e-50))
    n_perm_simple = int(merged.pop("n_perm_simple", 1000))
    if sample_size < 3:
        warnings.warn("sample_size is too small, so sample_size=3 is used.", RuntimeWarning, stacklevel=2)
        sample_size = 3
    if sample_size % 2 == 0:
        sample_size += 1
    eps = min(max(eps, 0.0), 1.0)

    simple_params = {
        **merged,
        "pathways": pathways,
        "stats": stats,
        "n_perm": n_perm_simple,
    }
    result = fgsea_simple(params=simple_params)
    if result.empty:
        return _empty_gsea_result(simple=False)
    result = result.drop(columns=["n_more_extreme"])
    result.insert(3, "log2err", multilevel_error(result["pval"].to_numpy(), sample_size=sample_size))
    hit_eps = result["pval"] < eps
    if hit_eps.any():
        result.loc[hit_eps, "pval"] = eps
        result.loc[hit_eps, "log2err"] = np.nan
        result["padj"] = p_adjust_bh(result["pval"])
        warnings.warn(
            f"For some pathways, in reality P-values are less than {eps}. Set eps=0 for better estimation.",
            RuntimeWarning,
            stacklevel=2,
        )
    return result[["pathway", "pval", "padj", "log2err", "es", "nes", "size", "leading_edge"]]


def fgsea(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run preranked GSEA.

    When ``n_perm`` is supplied, the simple permutation estimator is used.
    Otherwise the multilevel-compatible API is used.
    """

    merged = merge_params(params, allowed_keys=FGSEA_KEYS, context="fgsea", **kwargs)
    if pathways is not None:
        merged["pathways"] = pathways
    if stats is not None:
        merged["stats"] = stats
    if "n_perm" in merged and merged["n_perm"] is not None:
        simple_params = {key: value for key, value in merged.items() if key in FGSEA_SIMPLE_KEYS}
        return fgsea_simple(params=simple_params)
    multilevel_params = {key: value for key, value in merged.items() if key in FGSEA_MULTILEVEL_KEYS}
    return fgsea_multilevel(params=multilevel_params)


FGSEA_LABEL_KEYS = {
    "pathways",
    "expression",
    "labels",
    "n_perm",
    "min_size",
    "max_size",
    "n_jobs",
    "gsea_param",
    "score_type",
    "seed",
    "engine",
}


def fgsea_label(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    expression: Any = None,
    labels: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run label-permutation-style GSEA using gene-label correlations as ranks."""

    from ._utils import coerce_expression

    merged = merge_params(params, allowed_keys=FGSEA_LABEL_KEYS, context="fgsea_label", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "fgsea_label")
    expression = pop_required(merged, "expression", expression, "fgsea_label")
    labels = pop_required(merged, "labels", labels, "fgsea_label")
    frame = coerce_expression(expression)
    labels_array = pd.Series(labels).astype("category").cat.codes.to_numpy(dtype=float)
    if len(labels_array) != frame.shape[1]:
        raise ValueError("labels length must match the number of expression columns.")
    centered_labels = labels_array - labels_array.mean()
    expr = frame.to_numpy(dtype=float)
    expr = expr - expr.mean(axis=1, keepdims=True)
    denom = np.linalg.norm(expr, axis=1) * np.linalg.norm(centered_labels)
    stats = np.divide(expr @ centered_labels, denom, out=np.zeros(frame.shape[0], dtype=float), where=denom != 0)
    ranked = pd.Series(stats, index=frame.index)
    return fgsea_simple(pathways, ranked, params=merged)


COLLAPSE_KEYS = {"fgsea_res", "pathways", "stats", "pval_threshold", "overlap_threshold", "gsea_param"}


def collapse_pathways(
    fgsea_res: pd.DataFrame | None = None,
    pathways: Mapping[str, Iterable[Any]] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Collapse highly overlapping enriched pathways to representative pathways."""

    merged = merge_params(params, allowed_keys=COLLAPSE_KEYS, context="collapse_pathways", **kwargs)
    fgsea_res = pop_required(merged, "fgsea_res", fgsea_res, "collapse_pathways")
    pathways = pop_required(merged, "pathways", pathways, "collapse_pathways")
    _ = pop_required(merged, "stats", stats, "collapse_pathways")
    pval_threshold = float(merged.pop("pval_threshold", 0.05))
    overlap_threshold = float(merged.pop("overlap_threshold", 0.5))
    pathway_map = {name: set(genes) for name, genes in pathways.items()}
    ordered = fgsea_res.sort_values("pval", kind="mergesort")
    ordered = ordered[ordered["pval"] <= pval_threshold] if "pval" in ordered else ordered
    parent: dict[str, str | None] = {str(name): None for name in ordered["pathway"]}
    main: list[str] = []
    for name in ordered["pathway"].astype(str):
        if parent.get(name) is not None:
            continue
        main.append(name)
        genes = pathway_map.get(name, set())
        if not genes:
            continue
        for other in ordered["pathway"].astype(str):
            if other == name or parent.get(other) is not None:
                continue
            other_genes = pathway_map.get(other, set())
            denom = min(len(genes), len(other_genes)) or 1
            if len(genes & other_genes) / denom >= overlap_threshold:
                parent[other] = name
    return {"main_pathways": main, "parent_pathways": parent}
