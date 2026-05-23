# pyfgsea

`pyfgsea` is a Pythonic refactor of the `fgsea` fork in this repository. It
implements preranked GSEA, over-representation analysis, GESECA-style
co-regulation analysis, pathway I/O, and Matplotlib plotting helpers with an API
modeled after `PyOncoplot`.

The original R/C++ fork is preserved under
`python_refactor_goal_sources/fgsea_r/` for parity work and attribution.

## Install for Development

```bash
python3 -m pip install -e ".[test]"
```

## Quick Start

```python
import pandas as pd
from pyfgsea import fgsea, gmt_pathways, plot_enrichment

ranks = pd.Series(
    {"G1": 3.2, "G2": 2.7, "G3": 1.4, "G4": -0.5, "G5": -2.1, "G6": -2.8}
)
pathways = {
    "up": ["G1", "G2", "G3"],
    "down": ["G4", "G5", "G6"],
}

result = fgsea(pathways, ranks, min_size=2, max_size=5, seed=1)
print(result)

plot = plot_enrichment(pathways["up"], ranks, title="up")
plot.save("enrichment.png", dpi=120)
```

Every public API accepts reusable parameter dictionaries. Explicit keyword
arguments override values from `params`:

```python
params = {
    "pathways": pathways,
    "stats": ranks,
    "min_size": 2,
    "max_size": 5,
    "n_perm_simple": 500,
    "seed": 1,
}

result = fgsea(params=params, n_perm_simple=1000)
```

Unknown parameter names are rejected early with the allowed key list.

## Pythonic API

The public API uses snake_case names rather than R-style names:

| R `fgsea` name | Python `pyfgsea` name |
| --- | --- |
| `minSize` | `min_size` |
| `maxSize` | `max_size` |
| `scoreType` | `score_type` |
| `gseaParam` | `gsea_param` |
| `nperm` | `n_perm` |
| `nproc` | `n_jobs` |
| `sampleSize` | `sample_size` |
| `nPermSimple` | `n_perm_simple` |
| `leadingEdge` | `leading_edge` |

## Examples

The example gallery is driven by `examples/config.yaml`:

```bash
python3 examples/run_examples.py
```

Generated tables, plots, and GMT files are written to `examples/generated/`.
Run one preset with:

```bash
python3 examples/run_examples.py --preset quickstart_fgsea
```

## Engine Interface

Heavy APIs accept `engine="auto" | "python" | "compiled"`. The current v1
implementation ships the pure Python engine. `engine="auto"` resolves to Python
today, while `engine="compiled"` is reserved for the optional acceleration layer.

## Documentation

- [API reference](docs/api-reference.md)
- [Examples](docs/examples.md)
- [Migration from fgsea](docs/migration-from-fgsea.md)
