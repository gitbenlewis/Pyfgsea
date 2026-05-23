# Examples

The examples are configured in `examples/config.yaml` and run through
`examples/run_examples.py`.

```bash
python3 examples/run_examples.py
```

The config has five main sections:

- `output_dirs`: where generated tables, plots, and GMT files are written.
- `input_files`: deterministic ranks, pathways, expression, and gene lists.
- `default_params`: defaults merged into every run.
- `shared`: YAML anchors for reusable parameter blocks.
- `runs`: named presets with a renderer, output name, and params.

Generated artifacts are intentionally written under `examples/generated/` so the
source fixtures remain untouched.
