import numpy as np
import pandas as pd
import pytest

from pyfgsea import calc_gsea_stat, collapse_pathways, fgsea, fgsea_multilevel, fgsea_simple, multilevel_error


def test_calc_gsea_stat_numeric_cases():
    stats = pd.Series(list(range(10, -11, -1)), index=[f"G{i}" for i in range(21)])

    assert calc_gsea_stat(stats, list(range(5))) == pytest.approx(1)
    assert calc_gsea_stat(stats, list(range(15, 21))) == pytest.approx(-1)
    assert calc_gsea_stat(stats, [1, 3, 5, 7, 9]) == pytest.approx(0.71)


def test_calc_gsea_stat_leading_edge_with_gene_names():
    stats = pd.Series(list(range(10, -11, -1)), index=[f"G{i}" for i in range(21)])

    result = calc_gsea_stat(
        params={
            "stats": stats,
            "selected_stats": ["G10", "G0", "G1", "G2", "G3", "G4"],
            "return_leading_edge": True,
            "return_all_extremes": True,
        }
    )

    assert result["res"] == pytest.approx(calc_gsea_stat(stats, ["G10", "G0", "G1", "G2", "G3", "G4"]))
    assert "G0" in result["leading_edge"]
    assert "G10" not in result["leading_edge"]


def test_fgsea_validates_gene_names():
    ranks = pd.Series([1.0, 2.0, 3.0], index=["A", "A", "B"])
    with pytest.raises(ValueError, match="Duplicate"):
        fgsea_simple({"p": ["A", "B"]}, ranks, n_perm=10)

    ranks = pd.Series([1.0, np.inf], index=["A", "B"])
    with pytest.raises(ValueError, match="finite"):
        fgsea_simple({"p": ["A", "B"]}, ranks, n_perm=10)


def test_fgsea_handles_duplicate_pathway_genes_and_zero_pathways():
    ranks = pd.Series({"A": 3.0, "B": 2.0, "C": 1.0, "D": -1.0, "E": -2.0})

    result = fgsea_simple({"p1": ["A", "B"], "p2": ["A", "B", "A", "B"]}, ranks, n_perm=20, seed=1)
    assert result.set_index("pathway").loc["p1", "size"] == result.set_index("pathway").loc["p2", "size"]

    empty = fgsea_simple({"p": ["A", "B"]}, ranks, min_size=10, max_size=2, n_perm=10)
    assert empty.empty
    assert list(empty.columns) == ["pathway", "pval", "padj", "es", "nes", "n_more_extreme", "size", "leading_edge"]


def test_fgsea_params_kwargs_precedence_and_unknown_keys():
    ranks = pd.Series({"A": 3.0, "B": 2.0, "C": -1.0, "D": -2.0})
    pathways = {"p": ["A", "B"]}

    result = fgsea(params={"pathways": pathways, "stats": ranks, "n_perm_simple": 10, "seed": 1}, n_perm_simple=15)
    assert not result.empty

    with pytest.raises(ValueError, match="Unknown fgsea parameter"):
        fgsea(params={"pathways": pathways, "stats": ranks, "bad_name": True})


def test_fgsea_multilevel_and_error():
    ranks = pd.Series({"A": 3.0, "B": 2.0, "C": -1.0, "D": -2.0})
    pathways = {"p": ["A", "B"]}

    result = fgsea_multilevel(pathways, ranks, n_perm_simple=20, seed=1, eps=0)
    assert list(result.columns) == ["pathway", "pval", "padj", "log2err", "es", "nes", "size", "leading_edge"]
    assert multilevel_error(1e-4, sample_size=101) > 0


def test_collapse_pathways_overlap():
    ranks = pd.Series({"A": 3.0, "B": 2.0, "C": 1.0, "D": -1.0, "E": -2.0})
    pathways = {"p1": ["A", "B", "C"], "p2": ["A", "B", "C", "D"], "p3": ["D", "E"]}
    res = pd.DataFrame({"pathway": ["p1", "p2", "p3"], "pval": [0.001, 0.002, 0.5]})

    collapsed = collapse_pathways(res, pathways, ranks, pval_threshold=0.01, overlap_threshold=0.75)
    assert collapsed["main_pathways"][0] == "p1"
    assert collapsed["parent_pathways"]["p2"] == "p1"
