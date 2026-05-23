# API Reference

## GSEA

- `fgsea(pathways=None, stats=None, *, params=None, **kwargs)`
- `fgsea_simple(pathways=None, stats=None, *, params=None, **kwargs)`
- `fgsea_multilevel(pathways=None, stats=None, *, params=None, **kwargs)`
- `fgsea_label(pathways=None, expression=None, labels=None, *, params=None, **kwargs)`
- `calc_gsea_stat(stats=None, selected_stats=None, *, params=None, **kwargs)`
- `multilevel_error(pval=None, *, params=None, **kwargs)`
- `collapse_pathways(fgsea_res=None, pathways=None, stats=None, *, params=None, **kwargs)`

## ORA

- `fora(pathways=None, genes=None, universe=None, *, params=None, **kwargs)`
- `collapse_pathways_ora(fora_res=None, pathways=None, genes=None, universe=None, *, params=None, **kwargs)`

## GESECA

- `geseca(pathways=None, expression=None, *, params=None, **kwargs)`
- `geseca_simple(pathways=None, expression=None, *, params=None, **kwargs)`
- `collapse_pathways_geseca(geseca_res=None, pathways=None, expression=None, *, params=None, **kwargs)`

## Pathways

- `gmt_pathways(gmt_file=None, *, params=None, **kwargs)`
- `write_gmt_pathways(pathways=None, gmt_file=None, *, params=None, **kwargs)`
- `reactome_pathways(genes=None, *, params=None, **kwargs)`
- `map_ids_list(x=None, keys=None, column=None, keytype=None, *, params=None, **kwargs)`

## Plots

Plotting functions return `PlotResult(figure, axes, data)` and support
`result.save(path, **kwargs)`.

- `plot_enrichment_data(pathway=None, stats=None, *, params=None, **kwargs)`
- `plot_enrichment(pathway=None, stats=None, *, params=None, **kwargs)`
- `plot_gsea_table(pathways=None, stats=None, fgsea_res=None, *, params=None, **kwargs)`
- `plot_coregulation_profile(pathway=None, expression=None, *, params=None, **kwargs)`
- `plot_geseca_table(geseca_res=None, pathways=None, expression=None, *, params=None, **kwargs)`
- `plot_coregulation_profile_spatial(pathway=None, *, params=None, **kwargs)`
- `plot_coregulation_profile_reduction(pathway=None, *, params=None, **kwargs)`
- `plot_coregulation_profile_image(pathway=None, *, params=None, **kwargs)`
