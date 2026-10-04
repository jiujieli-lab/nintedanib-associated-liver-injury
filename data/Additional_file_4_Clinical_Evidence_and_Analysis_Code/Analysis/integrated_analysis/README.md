# Integrated revision analysis

The finalized revised inference tables are in `primary/`. The primary score uses the new donor-balanced phenotype-diffusion model, not numerically zero virtual-knockout ranks.

`uncorrected_vko_diagnostic/` reconstructs the earlier percentile-based workflow solely to document its numerical failure. Its target ranking must not be used as the revised result.

The root virtual-knockout audit and clinical-projection tables provide the data used for Figure 4 numeric validation. Filter the projection table to `numerically_nonzero == True` before biological interpretation.

Run commands (from a workspace with the input package and network outputs):

```bash
python audit_virtual_knockout_numerics.py --input-root INPUT_PACKAGE
python run_primary_integration.py --input-root INPUT_PACKAGE --network-dir NETWORK_DIRECTORY
python build_candidate_bridge_map.py --input-root INPUT_PACKAGE --network-dir NETWORK_DIRECTORY
python validate_primary_integration.py --input-root INPUT_PACKAGE --network-dir NETWORK_DIRECTORY
```

The independent GSE125975 liver dataset is intentionally excluded from algorithm fitting, feature selection, and weight selection.
