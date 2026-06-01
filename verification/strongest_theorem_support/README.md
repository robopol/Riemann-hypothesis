# Strongest Theorem Numerical Support Files

This folder contains reproducibility files for the numerical checks around
the strongest-beta and structured-candidate tests.  It documents only the
numerical audits and the data they generate.

## Contents

- `scripts/run_sigma_max_beta_test.py`
  Computes the last-prime beta envelope test and separates Euler's true
  `e^gamma` from the larger legacy GUI constant.

- `scripts/run_hcn_sa_bridge_test.py`
  Runs the HCN/SA-like structured log-profile search and tracks `log(N)`,
  `p_k`, `sigma(N)/N`, the Guy Robin index, and the bridge margins.

- `source_copies/`
  Contains preserved source scripts used as provenance for the headless
  numerical adaptations.

- `data/`
  Contains generated JSON outputs from the numerical audits.

- `figures/`
  Contains static charts generated from the audit outputs.

## Re-run

From this folder:

```powershell
py scripts/run_sigma_max_beta_test.py --max-count 1000000 --checkpoints 24
py scripts/run_hcn_sa_bridge_test.py --max-count 40000 --checkpoints 10
```

The scripts write their default outputs into `data/`.

## Status

These files are numerical support material.  They check scale, bookkeeping,
and candidate-profile behaviour; they are not a substitute for an analytic
proof.
