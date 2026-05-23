"""Run config-driven pyfgsea examples."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from pyfgsea import (  # noqa: E402
    fgsea,
    fgsea_multilevel,
    fgsea_simple,
    fora,
    geseca,
    gmt_pathways,
    plot_coregulation_profile,
    plot_enrichment,
    plot_geseca_table,
    plot_gsea_table,
    write_gmt_pathways,
)


CONFIG_PATH = Path(__file__).resolve().parent / "config.yaml"


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _repo_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def load_inputs(config: Mapping[str, Any]) -> dict[str, Any]:
    files = config["input_files"]
    ranks = pd.read_csv(_repo_path(files["ranks"]), sep="\t").set_index("gene")["stat"]
    pathways = gmt_pathways(_repo_path(files["pathways"]))
    genes = [line.strip() for line in _repo_path(files["genes"]).read_text(encoding="utf-8").splitlines() if line.strip()]
    expression = pd.read_csv(_repo_path(files["expression"]), sep="\t").set_index("gene")
    return {"ranks": ranks, "pathways": pathways, "genes": genes, "expression": expression}


def _table_out(config: Mapping[str, Any], name: str) -> Path:
    path = _repo_path(config["output_dirs"]["tables"]) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _plot_out(config: Mapping[str, Any], name: str) -> Path:
    path = _repo_path(config["output_dirs"]["plots"]) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _pathway_out(config: Mapping[str, Any], name: str) -> Path:
    path = _repo_path(config["output_dirs"]["pathways"]) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def run_one(name: str, run: Mapping[str, Any], config: Mapping[str, Any], inputs: Mapping[str, Any]) -> Path:
    params = _deep_merge(config.get("default_params", {}), run.get("params", {}))
    params.pop("run", None)
    save = params.pop("save", {}) or {}
    renderer = run["renderer"]
    output_name = run["output_name"]

    if renderer == "fgsea":
        result = fgsea(params={**params, "pathways": inputs["pathways"], "stats": inputs["ranks"]})
        output = _table_out(config, output_name)
        result.to_csv(output, sep="\t", index=False)
    elif renderer == "fgsea_simple":
        result = fgsea_simple(params={**params, "pathways": inputs["pathways"], "stats": inputs["ranks"]})
        output = _table_out(config, output_name)
        result.to_csv(output, sep="\t", index=False)
    elif renderer == "fgsea_multilevel":
        result = fgsea_multilevel(params={**params, "pathways": inputs["pathways"], "stats": inputs["ranks"]})
        output = _table_out(config, output_name)
        result.to_csv(output, sep="\t", index=False)
    elif renderer == "plot_enrichment":
        pathway = inputs["pathways"][params.pop("pathway")]
        params.pop("seed", None)
        result = plot_enrichment(params={**params, "pathway": pathway, "stats": inputs["ranks"]})
        output = _plot_out(config, output_name)
        result.save(output, dpi=int(save.get("dpi", 120)))
    elif renderer == "plot_gsea_table":
        analysis_params = dict(params)
        plot_params = dict(params)
        for key in ("top_n", "figure_size"):
            analysis_params.pop(key, None)
        for key in ("seed", "min_size", "max_size", "n_perm_simple", "eps", "sample_size", "n_jobs", "engine"):
            plot_params.pop(key, None)
        fgsea_res = fgsea(params={**analysis_params, "pathways": inputs["pathways"], "stats": inputs["ranks"]})
        result = plot_gsea_table(
            params={**plot_params, "pathways": inputs["pathways"], "stats": inputs["ranks"], "fgsea_res": fgsea_res}
        )
        output = _plot_out(config, output_name)
        result.save(output, dpi=int(save.get("dpi", 120)))
    elif renderer == "fora":
        params.pop("seed", None)
        result = fora(params={**params, "pathways": inputs["pathways"], "genes": inputs["genes"], "universe": inputs["ranks"].index})
        output = _table_out(config, output_name)
        result.to_csv(output, sep="\t", index=False)
    elif renderer == "geseca":
        result = geseca(params={**params, "pathways": inputs["pathways"], "expression": inputs["expression"]})
        output = _table_out(config, output_name)
        result.to_csv(output, sep="\t", index=False)
    elif renderer == "plot_coregulation_profile":
        pathway = inputs["pathways"][params.pop("pathway")]
        params.pop("seed", None)
        result = plot_coregulation_profile(params={**params, "pathway": pathway, "expression": inputs["expression"]})
        output = _plot_out(config, output_name)
        result.save(output, dpi=int(save.get("dpi", 120)))
    elif renderer == "plot_geseca_table":
        analysis_params = dict(params)
        plot_params = dict(params)
        for key in ("top_n", "figure_size", "min_limit", "max_limit", "titles"):
            analysis_params.pop(key, None)
        for key in ("seed", "min_size", "max_size", "n_perm_simple", "eps", "sample_size", "n_jobs", "engine"):
            plot_params.pop(key, None)
        geseca_res = geseca(params={**analysis_params, "pathways": inputs["pathways"], "expression": inputs["expression"]})
        result = plot_geseca_table(
            params={**plot_params, "geseca_res": geseca_res, "pathways": inputs["pathways"], "expression": inputs["expression"]}
        )
        output = _plot_out(config, output_name)
        result.save(output, dpi=int(save.get("dpi", 120)))
    elif renderer == "gmt_roundtrip":
        output = _pathway_out(config, output_name)
        write_gmt_pathways(inputs["pathways"], output)
    else:
        raise ValueError(f"Unknown renderer for {name}: {renderer}.")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preset", action="append", help="Run only the named preset. Can be supplied more than once.")
    args = parser.parse_args()
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
    inputs = load_inputs(config)
    selected = set(args.preset or [])
    for name, run in config["runs"].items():
        if selected and name not in selected:
            continue
        if not run.get("run", True):
            continue
        output = run_one(name, run, config, inputs)
        print(f"{name}: {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
