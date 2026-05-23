import numpy as np
import pandas as pd
import pytest

from pyfgsea import collapse_pathways_geseca, geseca, geseca_simple


def expression_frame():
    return pd.DataFrame(
        {
            "S1": [2.0, 1.8, 0.5, 0.4, 1.0],
            "S2": [2.2, 2.0, 0.4, 0.3, 0.8],
            "S3": [2.4, 2.1, 0.3, 0.2, 0.6],
            "S4": [0.6, 0.7, 1.8, 2.0, 1.1],
            "S5": [0.5, 0.6, 2.0, 2.2, 1.2],
        },
        index=["G1", "G2", "G3", "G4", "G5"],
    )


def test_geseca_reproducible():
    expression = expression_frame()
    pathways = {"up": ["G1", "G2"], "down": ["G3", "G4"]}

    one = geseca(pathways, expression, n_perm_simple=30, seed=1)
    two = geseca(pathways, expression, n_perm_simple=30, seed=1)

    assert one["pval"].tolist() == two["pval"].tolist()


def test_geseca_zero_pathways_and_columns():
    expression = expression_frame()
    empty = geseca({"p": ["G1", "G2"]}, expression, min_size=10, max_size=2, n_perm_simple=10, seed=1)

    assert empty.empty
    assert list(empty.columns) == ["pathway", "pct_var", "pval", "padj", "log2err", "size"]


def test_geseca_checks_gene_names_and_scale():
    expression = expression_frame()
    expression.index = ["G1", "G1", "G3", "G4", "G5"]

    with pytest.raises(ValueError, match="Duplicate"):
        geseca_simple({"p": ["G1", "G3"]}, expression, n_perm=10)

    expression = expression_frame()
    expression.loc["G5"] = 1.0
    with pytest.raises(ValueError, match="constant"):
        geseca_simple({"p": ["G1", "G2"]}, expression, scale=True, n_perm=10)


def test_geseca_eps_and_small_sample_warning():
    expression = expression_frame()
    pathways = {"up": ["G1", "G2"], "down": ["G3", "G4"]}

    with pytest.warns(RuntimeWarning, match="sample_size"):
        result = geseca(pathways, expression, n_perm_simple=20, seed=1, sample_size=1, eps=0)

    assert set(result.columns) == {"pathway", "pct_var", "pval", "padj", "log2err", "size"}


def test_collapse_pathways_geseca_overlap():
    expression = expression_frame()
    pathways = {"p1": ["G1", "G2"], "p2": ["G1", "G2", "G3"], "p3": ["G4", "G5"]}
    res = pd.DataFrame({"pathway": ["p1", "p2", "p3"], "pval": [0.01, 0.02, 0.5]})

    collapsed = collapse_pathways_geseca(res, pathways, expression, overlap_threshold=0.75)
    assert collapsed["parent_pathways"]["p2"] == "p1"
