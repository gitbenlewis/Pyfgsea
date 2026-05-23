"""Gene set co-regulation analysis."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any
import warnings

import numpy as np
import pandas as pd

from ._gsea import multilevel_error
from ._params import merge_params, pop_required
from ._utils import coerce_expression, p_adjust_bh, prepare_pathways, rng_from_seed, safe_log2err_from_count, validate_engine


GESECA_COMMON_KEYS = {
    "pathways",
    "expression",
    "min_size",
    "max_size",
    "center",
    "scale",
    "seed",
    "n_jobs",
    "engine",
}
GESECA_SIMPLE_KEYS = GESECA_COMMON_KEYS | {"n_perm"}
GESECA_KEYS = GESECA_COMMON_KEYS | {"sample_size", "eps", "n_perm_simple"}


def _scale_rows(frame: pd.DataFrame, *, center: bool, scale: bool) -> pd.DataFrame:
    values = frame.to_numpy(dtype=float)
    if center:
        values = values - values.mean(axis=1, keepdims=True)
    if scale:
        sd = values.std(axis=1, ddof=1, keepdims=True)
        if np.any(sd == 0):
            raise ValueError("Cannot rescale constant/zero gene expression rows to unit variance.")
        values = values / sd
    return pd.DataFrame(values, index=frame.index, columns=frame.columns)


def _geseca_score(values: np.ndarray, selected: list[int]) -> float:
    if not selected:
        return 0.0
    profile = values[selected, :].sum(axis=0)
    return float(np.sum(profile**2))


def _empty_geseca(simple: bool) -> pd.DataFrame:
    columns = ["pathway", "pct_var", "pval", "padj"]
    if not simple:
        columns.append("log2err")
    else:
        columns.append("n_more_extreme")
    columns.append("size")
    return pd.DataFrame({column: pd.Series(dtype=object) for column in columns})


def geseca_simple(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    expression: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run simple permutation-based GESECA."""

    merged = merge_params(params, allowed_keys=GESECA_SIMPLE_KEYS, context="geseca_simple", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "geseca_simple")
    expression = pop_required(merged, "expression", expression, "geseca_simple")
    min_size = int(merged.pop("min_size", 1))
    frame = coerce_expression(expression)
    max_size = int(merged.pop("max_size", max(frame.shape[0] - 1, 0)))
    center = bool(merged.pop("center", True))
    scale = bool(merged.pop("scale", False))
    n_perm = int(merged.pop("n_perm", 1000))
    seed = merged.pop("seed", None)
    n_jobs = int(merged.pop("n_jobs", 0))
    engine = validate_engine(merged.pop("engine", "auto"))
    if engine != "python":
        raise AssertionError("validate_engine should normalize to python.")
    if n_jobs not in (0, 1):
        warnings.warn("n_jobs is accepted for API compatibility; the Python engine currently runs serially.", RuntimeWarning)
    if n_perm < 1:
        raise ValueError("n_perm must be at least 1.")

    frame = _scale_rows(frame, center=center, scale=scale)
    filtered, sizes, _ = prepare_pathways(pathways, list(frame.index), min_size=min_size, max_size=max_size)
    if not filtered:
        return _empty_geseca(simple=True)

    values = frame.to_numpy(dtype=float)
    total_var = float(np.sum(values**2))
    rng = rng_from_seed(seed)
    rows: list[dict[str, Any]] = []
    for name, selected in filtered.items():
        size = sizes[name]
        actual = _geseca_score(values, selected)
        n_more_extreme = 0
        for _ in range(n_perm):
            random_selected = rng.choice(values.shape[0], size=size, replace=False).tolist()
            n_more_extreme += int(_geseca_score(values, random_selected) >= actual)
        pval = (n_more_extreme + 1) / (n_perm + 1)
        rows.append(
            {
                "pathway": name,
                "pct_var": actual / size / total_var * 100 if total_var and size else 0.0,
                "pval": pval,
                "n_more_extreme": n_more_extreme,
                "size": size,
            }
        )
    result = pd.DataFrame(rows)
    result["padj"] = p_adjust_bh(result["pval"])
    return result[["pathway", "pct_var", "pval", "padj", "n_more_extreme", "size"]].sort_values(
        "pval", kind="mergesort"
    ).reset_index(drop=True)


def geseca(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    expression: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """Run GESECA with the multilevel-compatible v1 Python estimator."""

    merged = merge_params(params, allowed_keys=GESECA_KEYS, context="geseca", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "geseca")
    expression = pop_required(merged, "expression", expression, "geseca")
    sample_size = int(merged.pop("sample_size", 101))
    eps = float(merged.pop("eps", 1e-50))
    n_perm_simple = int(merged.pop("n_perm_simple", 1000))
    if sample_size < 3:
        warnings.warn("sample_size is too small, so sample_size=3 is used.", RuntimeWarning, stacklevel=2)
        sample_size = 3
    if sample_size % 2 == 0:
        sample_size += 1
    eps = min(max(eps, 0.0), 1.0)
    result = geseca_simple(pathways, expression, params={**merged, "n_perm": n_perm_simple})
    if result.empty:
        return _empty_geseca(simple=False)
    result = result.drop(columns=["n_more_extreme"])
    result.insert(4, "log2err", multilevel_error(result["pval"].to_numpy(), sample_size=sample_size))
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
    return result[["pathway", "pct_var", "pval", "padj", "log2err", "size"]]


COLLAPSE_GESECA_KEYS = {
    "geseca_res",
    "pathways",
    "expression",
    "center",
    "scale",
    "eps",
    "check_depth",
    "n_jobs",
    "overlap_threshold",
}


def collapse_pathways_geseca(
    geseca_res: pd.DataFrame | None = None,
    pathways: Mapping[str, Iterable[Any]] | None = None,
    expression: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Collapse redundant GESECA pathways using gene overlap."""

    merged = merge_params(params, allowed_keys=COLLAPSE_GESECA_KEYS, context="collapse_pathways_geseca", **kwargs)
    geseca_res = pop_required(merged, "geseca_res", geseca_res, "collapse_pathways_geseca")
    pathways = pop_required(merged, "pathways", pathways, "collapse_pathways_geseca")
    _ = pop_required(merged, "expression", expression, "collapse_pathways_geseca")
    overlap_threshold = float(merged.pop("overlap_threshold", 0.5))
    pathway_map = {str(name): set(str(g) for g in genes) for name, genes in pathways.items()}
    ordered = geseca_res.sort_values("pval", kind="mergesort")
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
