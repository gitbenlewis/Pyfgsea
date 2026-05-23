# Migration From R fgsea

`pyfgsea` keeps the analysis concepts from `fgsea`, but it uses Pythonic naming
and standard Python data structures.

| R | Python |
| --- | --- |
| named numeric vector | `pandas.Series` indexed by gene |
| list of character vectors | `dict[str, list[str]]` |
| `data.table` result | `pandas.DataFrame` |
| `minSize` | `min_size` |
| `maxSize` | `max_size` |
| `scoreType` | `score_type` |
| `gseaParam` | `gsea_param` |
| `nperm` | `n_perm` |
| `nproc` | `n_jobs` |
| `sampleSize` | `sample_size` |
| `nPermSimple` | `n_perm_simple` |
| `leadingEdge` | `leading_edge` |

Seurat-specific plotting APIs are represented by AnnData/Scanpy-style helpers:
pass an object with `.X`, `.var_names`, and `.obsm`, or pass an expression
DataFrame plus explicit coordinates.
