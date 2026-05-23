from pathlib import Path

import pytest
from scipy.stats import hypergeom

from pyfgsea import fora, gmt_pathways, map_ids_list, write_gmt_pathways


def test_fora_hypergeometric_values():
    pathways = {
        "test_pathway_1": ["1", "2", "10"],
        "test_pathway_2": ["1", "2", "3", "4"],
        "test_pathway_3": ["6", "7", "8"],
    }
    genes = [str(i) for i in range(1, 6)]
    universe = [str(i) for i in range(1, 11)]

    result = fora(pathways, genes, universe)
    expected = sorted(
        [
            hypergeom.sf(2 - 1, 10, 3, 5),
            hypergeom.sf(4 - 1, 10, 4, 5),
            hypergeom.sf(0 - 1, 10, 3, 5),
        ]
    )
    assert result["pval"].tolist() == pytest.approx(expected)


def test_fora_drops_pathways_without_universe_overlap():
    result = fora({"p1": ["11"], "p2": ["1", "2"]}, genes=["1"], universe=["1", "2", "3"])
    assert result["pathway"].tolist() == ["p2"]


def test_gmt_roundtrip(tmp_path: Path):
    pathways = {"p1": ["A", "B"], "p2": ["C"]}
    path = tmp_path / "out.gmt"

    write_gmt_pathways(params={"pathways": pathways, "gmt_file": path})
    assert gmt_pathways(path) == pathways


def test_map_ids_list_mapping():
    mapped = map_ids_list({"A": "GeneA", "B": "GeneB"}, {"p": ["A", "B", "C"]}, missing=None)
    assert mapped == {"p": ["GeneA", "GeneB", None]}
