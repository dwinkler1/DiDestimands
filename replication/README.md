# DiD Estimand Companion — Replication Package

Companion code for the interactive applet at **https://didestimands.app/**.

## How it works

The applet is the data generator. The **"Panel (CSV)"** buttons on its
Replication tab export the simulated panel for the current slider settings as
`did_lab_panel.csv` (columns `unit, period, pair, treat, post, D, y`; one
button uses Tab 1's matched no-growth design, the other adds Tab 2's baseline
gap *b* and common growth *g*). Each script below reads that file, estimates
the four specifications, and prints estimates with 95% CIs. A reference panel
for the default settings (seed 4) is bundled, so every script runs out of the
box and can be checked against the reference values in its header before you
swap in your own export.

## Files

| File | Language | Estimation toolkit |
|---|---|---|
| `did_lab_replication.py` | Python ≥3.9 | `pyfixest`. **Reproducible:** `uv run did_lab_replication.py` — versions are pinned in a PEP 723 block at the top of the file. (Or `pip install numpy pandas scipy pyfixest`.) |
| `did_lab_replication.R`  | R ≥4.0      | `fixest`. **Reproducible:** the file opens with `renv::use("fixest@0.11.2")` to pin the version (needs `install.packages("renv")` once). |
| `did_lab_replication.jl` | Julia ≥1.6  | `DataFrames`, `FixedEffectModels`, `GLFixedEffectModels`, `FixedEffects`, `Distributions`, `DelimitedFiles` |
| `did_lab_panel.csv`      | —           | reference panel (default settings, seed 4) |

## What each script estimates

Four specifications, each with unit and period fixed effects and a single
Treat×Post regressor:

1. **Levels OLS**, reported as a percent of the counterfactual mean
   (treated pre-period mean × control post/pre ratio); its CI uses the delta
   method, because the scale is estimated and strongly negatively correlated
   with the coefficient (a fixed-scale CI would be ~2.5× too wide)
2. **log(1+Y) OLS**
3. **Weighted log(1+Y) OLS**, weights = unit's pre-period mean outcome
   (zero-weight units dropped — numerically identical)
4. **PPML** (Poisson pseudo-ML, log link), reported as exp(β)−1

All standard errors are CRV1, clustered at the matched-pair level, with
t critical values on (pairs−1) degrees of freedom. Two-way (unit × pair)
clustering — the paper's convention — is algebraically identical: units are
nested in pairs, so the CGM unit and intersection components cancel exactly.
The levels delta-method SE uses fixest's exact small-sample correction — the
finite-cluster factor G/(G−1) and the residual degrees-of-freedom factor
(N−1)/(N−K) with the nested-FE K (the unit FE absorbed within the pair
clusters are excluded), both inside the square root — so it matches `feols`/
`fepois` CRV1.

## Verification (bundled panel: defaults, seed 4)

```
n = 32,000   sum(y) = 23,551,023
levels  -0.03861113   [-4.4694%, -3.2529%]   (delta method, fixest-exact CRV1)
log1y   +0.00760246   [-0.3962%, +1.9167%]
wlog    -0.03704146   [-4.3714%, -3.0369%]
ppml    -0.03859994   [-4.4669%, -3.2492%]
```

- The applet's CSV export at the default sliders is **byte-identical** to the
  bundled panel, so "download → run" reproduces these numbers.
- **Python, R, and Julia** print digit-for-digit identical results (levels,
  log1y, wlog, ppml, and every CI end). PPML drops the same 60 separated
  (all-zero-unit) observations in all three.
- log/wlog/ppml interval ends can differ in the 3rd decimal across package
  versions (small-sample SE conventions); point estimates and the delta-method
  levels CI match exactly.

## Notes

- True estimand values (typical-unit %, population-total %) are properties of
  the data-generating process and are not recoverable from the CSV; read them
  off the applet's Tab 1 cards (defaults: −0.5% and −3.8%).
