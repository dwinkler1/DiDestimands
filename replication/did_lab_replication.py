# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy==1.23.5",
#     "pandas==1.5.3",
#     "scipy==1.10.1",
#     "pyfixest==0.50.1",
# ]
# ///
# Reproducible run:  uv run did_lab_replication.py  (uv reads the PEP 723 block above and
# builds a pinned, isolated environment). Versions are the ones this package was validated
# against; point estimates are version-independent.
"""
DiD estimation script (Python), companion to Winkler et al. (2026):
  https://doi.org/10.1007/s11129-026-09312-2

Estimates the companion's four specifications on the exported panel, all with unit
and period fixed effects and one Treat x Post regressor (D):
  levels OLS     reported as % of the counterfactual mean
                 (treated pre-mean x control post/pre ratio); delta-method CI
  log(1+Y) OLS
  weighted log(1+Y) OLS   (pre-period mean outcome weights)
  PPML           reported as exp(b) - 1

Input: did_lab_panel.csv. Download it from the applet ("Panel (CSV)" buttons
on the Replication tab; any slider settings) or use the bundled copy
(default settings, seed 4). Reference values for the bundled panel:
  n = 32,000   sum(y) = 23,551,023
  levels -0.03861113 [-4.4694%, -3.2529%]   log(1+Y) +0.00760246 [-0.3962%, +1.9167%]
  w-log  -0.03704146 [-4.3714%, -3.0369%]   PPML     -0.03859994 [-4.4669%, -3.2492%]
  (levels CI uses the fixest-exact CRV1 correction; other CI ends can shift in the 3rd decimal across versions)
True estimand values are properties of the DGP, not recoverable from the CSV;
read them off the applet's Tab 1 cards.

Requires: numpy, pandas, scipy, pyfixest. Best run with `uv run` (PEP 723 block at top).
"""

import numpy as np
import pandas as pd
import pyfixest as pf
from pyfixest.estimation import demean
from scipy.stats import t as t_dist

df = pd.read_csv("did_lab_panel.csv")
print(f"rows = {len(df)}   sum(y) = {int(df.y.sum())}")

# scale for the levels estimate: treated pre-mean x control post/pre ratio
mean_treat_pre = df.loc[(df.treat == 1) & (df.post == 0), "y"].mean()
mean_ctrl_pre  = df.loc[(df.treat == 0) & (df.post == 0), "y"].mean()
mean_ctrl_post = df.loc[(df.treat == 0) & (df.post == 1), "y"].mean()
scale = mean_treat_pre * mean_ctrl_post / mean_ctrl_pre

n_pairs = df["pair"].nunique()
tcrit   = t_dist.ppf(0.975, df=n_pairs - 1)

df["log_y"] = np.log1p(df.y)
w_pre = df[df.post == 0].groupby("unit")["y"].mean().rename("w_pre")
df = df.merge(w_pre, on="unit")

# SEs: CRV1 clustered by matched pair. With 1:1 matching WITHOUT replacement
# (this panel), each unit belongs to exactly one pair, so two-way unit-and-pair
# clustering is identical to pair clustering. If you match WITH replacement
# (controls reused across pairs), units no longer nest in pairs, so cluster
# two-way instead: vcov={"CRV1": "unit+pair"}.
V = {"CRV1": "pair"}
m_levels = pf.feols("y ~ D | unit + period",     data=df, vcov=V)
m_log    = pf.feols("log_y ~ D | unit + period", data=df, vcov=V)
m_wlog   = pf.feols("log_y ~ D | unit + period", data=df[df.w_pre > 0],
                    weights="w_pre", vcov=V)  # pyfixest requires strictly
                    # positive weights, so drop zero-weight units (they carry
                    # no information either way)
m_ppml   = pf.fepois("y ~ D | unit + period",    data=df, vcov=V)
# (fepois may note that all-zero units were dropped for separation; this is expected)
b  = lambda m: m.coef()["D"]
se = lambda m: m.tidy()["Std. Error"]["D"]

# delta-method CI for the scaled levels effect theta = b/scale: the scale is
# estimated from the same sample and strongly negatively correlated with b,
# so a fixed-scale CI would be ~2.5x too wide.
D_dem, _ = demean(df[["D"]].to_numpy(float),
                  df[["unit", "period"]].to_numpy(), np.ones(len(df)))
D_dem = D_dem[:, 0]
denom = D_dem @ D_dem
resid_lv   = m_levels.resid()
score_pair = np.bincount(df["pair"], weights=D_dem * resid_lv,
                         minlength=n_pairs)
pair_mean    = lambda rows: df[rows].groupby("pair")["y"].mean().to_numpy()
pm_treat_pre = pair_mean((df.treat == 1) & (df.post == 0))
pm_ctrl_pre  = pair_mean((df.treat == 0) & (df.post == 0))
pm_ctrl_post = pair_mean((df.treat == 0) & (df.post == 1))
theta = b(m_levels) / scale
psi = (score_pair / (denom * scale)
       - (theta / n_pairs) * (
             (pm_treat_pre - pm_treat_pre.mean()) / pm_treat_pre.mean()
           - (pm_ctrl_pre  - pm_ctrl_pre.mean())  / pm_ctrl_pre.mean()
           + (pm_ctrl_post - pm_ctrl_post.mean()) / pm_ctrl_post.mean()))
# CRV1, fixest-exact: BOTH small-sample factors multiply the VARIANCE (inside the sqrt),
# not the SE. Nested K = period FE + D (unit FE absorbed within the pair clusters, excluded).
K_dof = 1 + df["period"].nunique()
ssc   = (n_pairs / (n_pairs - 1)) * ((len(df) - 1) / (len(df) - K_dof))
se_theta = np.sqrt((psi**2).sum() * ssc)

ci = lambda est, s: f"[{100*(est - tcrit*s):+.4f}%, {100*(est + tcrit*s):+.4f}%]"
print(f"levels  {theta:+.8f}   {ci(theta, se_theta)}   (delta method)")
print(f"log1y   {b(m_log):+.8f}   {ci(b(m_log), se(m_log))}")
print(f"wlog    {b(m_wlog):+.8f}   {ci(b(m_wlog), se(m_wlog))}")
d = b(m_ppml)
print(f"ppml    {np.exp(d)-1:+.8f}   [{100*(np.exp(d - tcrit*se(m_ppml))-1):+.4f}%, "
      f"{100*(np.exp(d + tcrit*se(m_ppml))-1):+.4f}%]")


# Standard regression table (raw coefficients).
# NOTE: levels is in outcome units here (not the scaled % above), PPML is the
# log coefficient (not exp(b)-1), and SEs use pyfixest's default convention.
print()
pf.etable([m_levels, m_log, m_wlog, m_ppml], type="md")
