"""Matplotlib plotting helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._gsea import calc_gsea_stat
from ._geseca import _scale_rows
from ._params import merge_params, pop_required
from ._result import PlotResult
from ._utils import coerce_expression, prepare_ranked_stats


PLOT_ENRICHMENT_DATA_KEYS = {"pathway", "stats", "gsea_param", "stats_gene_col", "stats_value_col"}
PLOT_ENRICHMENT_KEYS = PLOT_ENRICHMENT_DATA_KEYS | {"ticks_size", "title", "figure_size"}
PLOT_GSEA_TABLE_KEYS = {"pathways", "stats", "fgsea_res", "gsea_param", "top_n", "figure_size"}
PROFILE_KEYS = {"pathway", "expression", "center", "scale", "titles", "conditions", "title", "figure_size"}
GESECA_TABLE_KEYS = {
    "geseca_res",
    "pathways",
    "expression",
    "center",
    "scale",
    "titles",
    "top_n",
    "figure_size",
    "min_limit",
    "max_limit",
}
SPATIAL_KEYS = {
    "pathway",
    "adata",
    "object",
    "expression",
    "coordinates",
    "basis",
    "title",
    "figure_size",
    "min_limit",
    "max_limit",
}


def plot_enrichment_data(
    pathway: Iterable[Any] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """Return data used to draw a GSEA enrichment plot."""

    merged = merge_params(params, allowed_keys=PLOT_ENRICHMENT_DATA_KEYS, context="plot_enrichment_data", **kwargs)
    pathway = pop_required(merged, "pathway", pathway, "plot_enrichment_data")
    stats = pop_required(merged, "stats", stats, "plot_enrichment_data")
    gsea_param = float(merged.pop("gsea_param", 1.0))
    stats_gene_col = merged.pop("stats_gene_col", None)
    stats_value_col = merged.pop("stats_value_col", None)
    ranked = prepare_ranked_stats(stats, gsea_param=gsea_param, gene_col=stats_gene_col, value_col=stats_value_col)
    genes = list(ranked.index)
    positions = sorted({genes.index(str(gene)) for gene in pathway if str(gene) in ranked.index})
    stat = calc_gsea_stat(
        ranked,
        positions,
        gsea_param=gsea_param,
        return_all_extremes=True,
        params=None,
    )
    assert isinstance(stat, dict)
    values = ranked.to_numpy()
    signs = np.sign(values)
    adjusted = np.abs(values) ** gsea_param
    n = len(values)
    xs = np.ravel(np.column_stack((np.array(positions), np.array(positions) + 1))) if positions else np.array([])
    ys = np.ravel(np.column_stack((stat["bottoms"], stat["tops"]))) if positions else np.array([])
    curve = pd.DataFrame({"rank": np.concatenate(([0], xs, [n])), "es": np.concatenate(([0.0], ys, [0.0]))})
    ticks = pd.DataFrame({"rank": np.array(positions) + 1, "stat": signs[positions] * adjusted[positions]})
    stats_frame = pd.DataFrame({"rank": np.arange(1, n + 1), "stat": signs * adjusted})
    return {
        "curve": curve,
        "ticks": ticks,
        "stats": stats_frame,
        "pos_es": float(np.max(stat["tops"])) if positions else 0.0,
        "neg_es": float(np.min(stat["bottoms"])) if positions else 0.0,
        "spread_es": float(np.max(stat["tops"]) - np.min(stat["bottoms"])) if positions else 0.0,
        "max_abs_stat": float(np.max(np.abs(adjusted))) if len(adjusted) else 0.0,
    }


def plot_enrichment(
    pathway: Iterable[Any] | None = None,
    stats: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Draw a GSEA enrichment plot."""

    merged = merge_params(params, allowed_keys=PLOT_ENRICHMENT_KEYS, context="plot_enrichment", **kwargs)
    ticks_size = float(merged.pop("ticks_size", 0.5))
    title = merged.pop("title", None)
    figure_size = tuple(merged.pop("figure_size", (7, 4)))
    data = plot_enrichment_data(pathway, stats, params=merged)
    fig, ax = plt.subplots(figsize=figure_size)
    ax.plot(data["curve"]["rank"], data["curve"]["es"], color="#1B9E77", linewidth=1.8)
    spread = data["spread_es"] or 1.0
    ax.vlines(data["ticks"]["rank"], -spread / 16, spread / 16, color="black", linewidth=ticks_size)
    ax.axhline(data["pos_es"], color="#D95F02", linestyle="--", linewidth=1)
    ax.axhline(data["neg_es"], color="#D95F02", linestyle="--", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("rank")
    ax.set_ylabel("enrichment score")
    if title:
        ax.set_title(str(title))
    fig.tight_layout()
    return PlotResult(fig, ax, data)


def plot_gsea_table(
    pathways: Mapping[str, Iterable[Any]] | None = None,
    stats: Any = None,
    fgsea_res: pd.DataFrame | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Draw a compact table of ranked-gene barcode plots and GSEA values."""

    merged = merge_params(params, allowed_keys=PLOT_GSEA_TABLE_KEYS, context="plot_gsea_table", **kwargs)
    pathways = pop_required(merged, "pathways", pathways, "plot_gsea_table")
    stats = pop_required(merged, "stats", stats, "plot_gsea_table")
    fgsea_res = pop_required(merged, "fgsea_res", fgsea_res, "plot_gsea_table")
    gsea_param = float(merged.pop("gsea_param", 1.0))
    top_n = int(merged.pop("top_n", min(15, len(fgsea_res))))
    figure_size = tuple(merged.pop("figure_size", (9, max(3, top_n * 0.45 + 1))))
    ranked = prepare_ranked_stats(stats, gsea_param=gsea_param)
    order = fgsea_res.head(top_n)["pathway"].astype(str).tolist()
    fig, axes = plt.subplots(len(order), 1, figsize=figure_size, sharex=True)
    if len(order) == 1:
        axes = np.array([axes])
    for ax, name in zip(axes, order):
        genes = pathways.get(name, [])
        positions = [i for i, gene in enumerate(ranked.index) if gene in set(str(g) for g in genes)]
        ax.vlines(positions, 0, 1, color="#333333", linewidth=0.6)
        row = fgsea_res[fgsea_res["pathway"].astype(str) == name].iloc[0]
        ax.set_ylabel(str(name), rotation=0, ha="right", va="center")
        ax.text(1.01, 0.5, f"NES {row.get('nes', np.nan):.2f}  p {row.get('pval', np.nan):.2g}", transform=ax.transAxes)
        ax.set_yticks([])
    axes[-1].set_xlabel("ranked genes")
    fig.tight_layout()
    return PlotResult(fig, axes, {"pathways": order})


def plot_coregulation_profile(
    pathway: Iterable[Any] | None = None,
    expression: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Plot expression profiles for a gene set."""

    merged = merge_params(params, allowed_keys=PROFILE_KEYS, context="plot_coregulation_profile", **kwargs)
    pathway = pop_required(merged, "pathway", pathway, "plot_coregulation_profile")
    expression = pop_required(merged, "expression", expression, "plot_coregulation_profile")
    center = bool(merged.pop("center", True))
    scale = bool(merged.pop("scale", False))
    titles = merged.pop("titles", None)
    conditions = merged.pop("conditions", None)
    title = merged.pop("title", None)
    figure_size = tuple(merged.pop("figure_size", (9, 4)))
    frame = _scale_rows(coerce_expression(expression), center=center, scale=scale)
    genes = [str(gene) for gene in pathway if str(gene) in frame.index]
    selected = frame.loc[genes]
    x = np.arange(selected.shape[1])
    labels = [str(x) for x in (titles if titles is not None else selected.columns)]
    fig, ax = plt.subplots(figsize=figure_size)
    for _, row in selected.iterrows():
        ax.plot(x, row.to_numpy(dtype=float), color="#6C8EBF", alpha=0.25, linewidth=0.8)
        ax.scatter(x, row.to_numpy(dtype=float), color="#6C8EBF", alpha=0.25, s=10)
    profile = selected.mean(axis=0).to_numpy(dtype=float) if len(selected) else np.zeros(len(labels))
    if conditions is not None:
        codes = pd.Series(conditions).astype("category").cat.codes.to_numpy()
        ax.scatter(x, profile, c=codes, cmap="tab10", edgecolor="black", s=36, zorder=3)
    else:
        ax.scatter(x, profile, color="black", s=36, zorder=3)
    ax.plot(x, profile, color="#13242A", linewidth=1.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("expression")
    if title:
        ax.set_title(str(title))
    fig.tight_layout()
    return PlotResult(fig, ax, {"genes": genes, "profile": profile})


def plot_geseca_table(
    geseca_res: pd.DataFrame | None = None,
    pathways: Mapping[str, Iterable[Any]] | None = None,
    expression: Any = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Draw a compact GESECA pathway profile heatmap."""

    merged = merge_params(params, allowed_keys=GESECA_TABLE_KEYS, context="plot_geseca_table", **kwargs)
    geseca_res = pop_required(merged, "geseca_res", geseca_res, "plot_geseca_table")
    pathways = pop_required(merged, "pathways", pathways, "plot_geseca_table")
    expression = pop_required(merged, "expression", expression, "plot_geseca_table")
    center = bool(merged.pop("center", True))
    scale = bool(merged.pop("scale", False))
    titles = merged.pop("titles", None)
    top_n = int(merged.pop("top_n", min(10, len(geseca_res))))
    figure_size = tuple(merged.pop("figure_size", (9, max(3, top_n * 0.45 + 1))))
    min_limit = float(merged.pop("min_limit", -3))
    max_limit = float(merged.pop("max_limit", 3))
    frame = _scale_rows(coerce_expression(expression), center=center, scale=scale)
    order = geseca_res.head(top_n)["pathway"].astype(str).tolist()
    profiles = []
    for name in order:
        genes = [str(g) for g in pathways.get(name, []) if str(g) in frame.index]
        profile = frame.loc[genes].sum(axis=0).to_numpy(dtype=float) if genes else np.zeros(frame.shape[1])
        sd = profile.std(ddof=1)
        profiles.append((profile - profile.mean()) / sd if sd else profile * 0)
    data = np.vstack(profiles) if profiles else np.zeros((0, frame.shape[1]))
    fig, ax = plt.subplots(figsize=figure_size)
    image = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=min_limit, vmax=max_limit)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xticks(np.arange(frame.shape[1]))
    ax.set_xticklabels([str(x) for x in (titles if titles is not None else frame.columns)], rotation=45, ha="right")
    fig.colorbar(image, ax=ax, label="z-score")
    fig.tight_layout()
    return PlotResult(fig, ax, pd.DataFrame(data, index=order, columns=frame.columns))


def _score_cells_from_adata(pathway: Iterable[Any], adata: Any) -> tuple[np.ndarray, list[str]]:
    genes = [str(gene) for gene in pathway]
    var_names = [str(gene) for gene in getattr(adata, "var_names")]
    positions = [var_names.index(gene) for gene in genes if gene in var_names]
    if not positions:
        raise ValueError("No pathway genes were found in adata.var_names.")
    x = adata.X
    if hasattr(x, "toarray"):
        x = x.toarray()
    values = np.asarray(x, dtype=float)[:, positions]
    return values.sum(axis=1) / np.sqrt(len(positions)), [var_names[i] for i in positions]


def _scatter_profile(
    pathway: Iterable[Any] | None,
    *,
    params: Mapping[str, Any] | None,
    context: str,
    default_basis: str,
    **kwargs: Any,
) -> PlotResult:
    merged = merge_params(params, allowed_keys=SPATIAL_KEYS, context=context, **kwargs)
    pathway = pop_required(merged, "pathway", pathway, context)
    adata = merged.pop("adata", None)
    adata = merged.pop("object", None) if adata is None else adata
    expression = merged.pop("expression", None)
    coordinates = merged.pop("coordinates", None)
    basis = merged.pop("basis", default_basis)
    title = merged.pop("title", None)
    figure_size = tuple(merged.pop("figure_size", (6, 5)))
    min_limit = merged.pop("min_limit", None)
    max_limit = merged.pop("max_limit", None)
    if adata is not None:
        scores, genes = _score_cells_from_adata(pathway, adata)
        if not hasattr(adata, "obsm") or basis not in adata.obsm:
            raise ValueError(f"adata.obsm must contain coordinates for basis {basis!r}.")
        coords = np.asarray(adata.obsm[basis])
    else:
        expression = pop_required(merged, "expression", expression, context)
        coordinates = pop_required(merged, "coordinates", coordinates, context)
        frame = coerce_expression(expression)
        genes = [str(gene) for gene in pathway if str(gene) in frame.index]
        scores = frame.loc[genes].sum(axis=0).to_numpy(dtype=float) / np.sqrt(len(genes))
        coords = np.asarray(coordinates, dtype=float)
    fig, ax = plt.subplots(figsize=figure_size)
    points = ax.scatter(coords[:, 0], coords[:, 1], c=scores, cmap="RdBu_r", vmin=min_limit, vmax=max_limit)
    fig.colorbar(points, ax=ax, label="pathway score")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if title:
        ax.set_title(str(title))
    fig.tight_layout()
    return PlotResult(fig, ax, {"scores": scores, "genes": genes})


def plot_coregulation_profile_spatial(
    pathway: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Plot pathway scores over AnnData ``obsm['spatial']`` or supplied coordinates."""

    return _scatter_profile(pathway, params=params, context="plot_coregulation_profile_spatial", default_basis="spatial", **kwargs)


def plot_coregulation_profile_reduction(
    pathway: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Plot pathway scores over AnnData reduction coordinates such as ``X_umap``."""

    return _scatter_profile(pathway, params=params, context="plot_coregulation_profile_reduction", default_basis="X_umap", **kwargs)


def plot_coregulation_profile_image(
    pathway: Iterable[Any] | None = None,
    *,
    params: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> PlotResult:
    """Alias for spatial image-style pathway score plotting."""

    return _scatter_profile(pathway, params=params, context="plot_coregulation_profile_image", default_basis="spatial", **kwargs)
