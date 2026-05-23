from pathlib import Path
import subprocess
import sys

import matplotlib
import pandas as pd

matplotlib.use("Agg")

from pyfgsea import (  # noqa: E402
    fgsea,
    geseca,
    plot_coregulation_profile,
    plot_enrichment,
    plot_enrichment_data,
    plot_geseca_table,
    plot_gsea_table,
)


def ranks():
    return pd.Series({"G1": 3.0, "G2": 2.0, "G3": 1.0, "G4": -1.0, "G5": -2.0, "G6": -3.0})


def pathways():
    return {"up": ["G1", "G2", "G3"], "down": ["G4", "G5", "G6"]}


def expression():
    return pd.DataFrame(
        {"S1": [2, 2, 1, 0, 0, 1], "S2": [3, 3, 1, 0, 0, 1], "S3": [0, 0, 1, 3, 3, 1]},
        index=["G1", "G2", "G3", "G4", "G5", "G6"],
    )


def test_plot_smoke(tmp_path: Path):
    fgsea_res = fgsea(pathways(), ranks(), n_perm_simple=20, seed=1)
    geseca_res = geseca(pathways(), expression(), n_perm_simple=20, seed=1)

    for result, name in [
        (plot_enrichment(pathways()["up"], ranks()), "enrichment.png"),
        (plot_gsea_table(pathways(), ranks(), fgsea_res), "gsea_table.png"),
        (plot_coregulation_profile(pathways()["up"], expression()), "profile.png"),
        (plot_geseca_table(geseca_res, pathways(), expression()), "geseca_table.png"),
    ]:
        output = tmp_path / name
        result.save(output)
        assert output.exists()

    data = plot_enrichment_data(pathways()["up"], ranks())
    assert {"curve", "ticks", "stats", "pos_es", "neg_es"}.issubset(data)


def test_config_driven_examples_selected():
    root = Path(__file__).resolve().parents[1]
    subprocess.run(
        [sys.executable, "examples/run_examples.py", "--preset", "quickstart_fgsea", "--preset", "gmt_roundtrip"],
        cwd=root,
        check=True,
    )
    assert (root / "examples/generated/tables/quickstart_fgsea.tsv").exists()
    assert (root / "examples/generated/pathways/roundtrip.gmt").exists()
