"""Pythonic fast gene set enrichment analysis."""

from ._geseca import collapse_pathways_geseca, geseca, geseca_simple
from ._gsea import (
    calc_gsea_stat,
    collapse_pathways,
    fgsea,
    fgsea_label,
    fgsea_multilevel,
    fgsea_simple,
    multilevel_error,
)
from ._ora import collapse_pathways_ora, fora
from ._pathways import gmt_pathways, map_ids_list, reactome_pathways, write_gmt_pathways
from ._plotting import (
    plot_coregulation_profile,
    plot_coregulation_profile_image,
    plot_coregulation_profile_reduction,
    plot_coregulation_profile_spatial,
    plot_enrichment,
    plot_enrichment_data,
    plot_geseca_table,
    plot_gsea_table,
)
from ._result import PlotResult

__all__ = [
    "PlotResult",
    "calc_gsea_stat",
    "collapse_pathways",
    "collapse_pathways_geseca",
    "collapse_pathways_ora",
    "fgsea",
    "fgsea_label",
    "fgsea_multilevel",
    "fgsea_simple",
    "fora",
    "geseca",
    "geseca_simple",
    "gmt_pathways",
    "map_ids_list",
    "multilevel_error",
    "plot_coregulation_profile",
    "plot_coregulation_profile_image",
    "plot_coregulation_profile_reduction",
    "plot_coregulation_profile_spatial",
    "plot_enrichment",
    "plot_enrichment_data",
    "plot_geseca_table",
    "plot_gsea_table",
    "reactome_pathways",
    "write_gmt_pathways",
]
