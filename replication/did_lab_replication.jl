# ==============================================================================
# DiD estimation script (Julia), companion to Winkler et al. (2026):
#   https://doi.org/10.1007/s11129-026-09312-2
#
# Estimates the companion's four specifications on the exported panel, all with unit
# and period fixed effects and one Treat x Post regressor (D):
#   levels OLS     reported as % of the counterfactual mean
#                  (treated pre-mean x control post/pre ratio); delta-method CI
#   log(1+Y) OLS
#   weighted log(1+Y) OLS   (pre-period mean outcome weights)
#   PPML           reported as exp(b) - 1
#
# Input: did_lab_panel.csv. Download it from the applet ("Panel (CSV)" buttons
# on the Replication tab; any slider settings) or use the bundled copy
# (default settings, seed 4). Reference values for the bundled panel:
#   n = 32,000   sum(y) = 23,551,023
#   levels -0.03861113 [-4.4694%, -3.2529%]   log(1+Y) +0.00760246 [-0.3962%, +1.9167%]
#   w-log  -0.03704146 [-4.3714%, -3.0369%]   PPML     -0.03859994 [-4.4669%, -3.2492%]
#   (levels CI uses the fixest-exact CRV1 correction; other CI ends can shift in the 3rd decimal across versions)
# True estimand values are properties of the DGP, not recoverable from the CSV;
# read them off the applet's Tab 1 cards.
#
# Requires: import Pkg; Pkg.add(["DelimitedFiles","DataFrames","Distributions",
#   "FixedEffects","FixedEffectModels","GLFixedEffectModels","RegressionTables"])
# ==============================================================================

using DelimitedFiles, DataFrames, Statistics, Printf
using FixedEffectModels, GLFixedEffectModels
using FixedEffects: FixedEffect, solve_residuals!
using RegressionTables
using Distributions: TDist

M, hdr = readdlm("did_lab_panel.csv", ',', Int; header = true)
df = DataFrame(M, vec(Symbol.(hdr)))
@printf("rows = %d   sum(y) = %d\n", nrow(df), sum(df.y))

# scale for the levels estimate: treated pre-mean x control post/pre ratio
mean_treat_pre = mean(df.y[(df.treat .== 1) .& (df.post .== 0)])
mean_ctrl_pre  = mean(df.y[(df.treat .== 0) .& (df.post .== 0)])
mean_ctrl_post = mean(df.y[(df.treat .== 0) .& (df.post .== 1)])
scale = mean_treat_pre * mean_ctrl_post / mean_ctrl_pre

n_pairs = length(unique(df.pair))
tcrit   = quantile(TDist(n_pairs - 1), 0.975)

df.log_y = log1p.(df.y)
w_pre = combine(groupby(df[df.post .== 0, :], :unit), :y => mean => :w_pre)
df = innerjoin(df, w_pre, on = :unit)

# SEs: CRV1 clustered by matched pair. With 1:1 matching WITHOUT replacement
# (this panel), each unit belongs to exactly one pair, so two-way unit-and-pair
# clustering is identical to pair clustering. If you match WITH replacement
# (controls reused across pairs), units no longer nest in pairs, so cluster
# two-way instead: Vcov.cluster(:unit, :pair).
m_levels = reg(df, @formula(y ~ D + fe(unit) + fe(period)), Vcov.cluster(:pair),
               save = :residuals)
m_log    = reg(df, @formula(log_y ~ D + fe(unit) + fe(period)), Vcov.cluster(:pair))
m_wlog   = reg(df[df.w_pre .> 0, :], @formula(log_y ~ D + fe(unit) + fe(period)),
               Vcov.cluster(:pair), weights = :w_pre)  # zero-weight units drop out anyway
m_ppml   = nlreg(df, @formula(y ~ D + fe(unit) + fe(period)),
                 Poisson(), LogLink(), Vcov.cluster(:pair),
                 separation = [:fe])   # drop all-zero units, like the other languages
bof(m) = coef(m)[findfirst(==("D"), coefnames(m))]
sof(m) = stderror(m)[findfirst(==("D"), coefnames(m))]

# delta-method CI for the scaled levels effect theta = b/scale: the scale is
# estimated from the same sample and strongly negatively correlated with b,
# so a fixed-scale CI would be ~2.5x too wide.
D_dem = solve_residuals!(Float64.(df.D),
                         [FixedEffect(df.unit), FixedEffect(df.period)])[1]
resid_lv = Float64.(residuals(m_levels))
denom = sum(abs2, D_dem)
df.score = D_dem .* resid_lv
score_pair = sort!(combine(groupby(df, :pair), :score => sum => :s), :pair).s
pair_mean(rows) = sort!(combine(groupby(df[rows, :], :pair), :y => mean => :m), :pair).m
pm_treat_pre = pair_mean((df.treat .== 1) .& (df.post .== 0))
pm_ctrl_pre  = pair_mean((df.treat .== 0) .& (df.post .== 0))
pm_ctrl_post = pair_mean((df.treat .== 0) .& (df.post .== 1))
theta = bof(m_levels) / scale
psi = score_pair ./ (denom * scale) .-
      (theta / n_pairs) .* ((pm_treat_pre .- mean(pm_treat_pre)) ./ mean(pm_treat_pre) .-
                            (pm_ctrl_pre  .- mean(pm_ctrl_pre))  ./ mean(pm_ctrl_pre)  .+
                            (pm_ctrl_post .- mean(pm_ctrl_post)) ./ mean(pm_ctrl_post))
# CRV1, fixest-exact: BOTH small-sample factors multiply the VARIANCE (inside the sqrt).
# Nested K = period FE + D (unit FE absorbed within the pair clusters, excluded).
K_dof    = 1 + length(unique(df.period))
ssc      = (n_pairs / (n_pairs - 1)) * ((nrow(df) - 1) / (nrow(df) - K_dof))
se_theta = sqrt(sum(abs2, psi) * ssc)

ci(est, s) = @sprintf("[%+.4f%%, %+.4f%%]", 100*(est - tcrit*s), 100*(est + tcrit*s))
@printf("levels  %+.8f   %s   (delta method)\n", theta, ci(theta, se_theta))
@printf("log1y   %+.8f   %s\n", bof(m_log),  ci(bof(m_log),  sof(m_log)))
@printf("wlog    %+.8f   %s\n", bof(m_wlog), ci(bof(m_wlog), sof(m_wlog)))
d = bof(m_ppml)
@printf("ppml    %+.8f   [%+.4f%%, %+.4f%%]\n", exp(d) - 1,
        100*(exp(d - tcrit*sof(m_ppml)) - 1), 100*(exp(d + tcrit*sof(m_ppml)) - 1))


# Standard regression table (raw coefficients).
# NOTE: levels is in outcome units here (not the scaled % above), PPML is the
# log coefficient (not exp(b)-1), and SEs use the package default convention.
println(regtable(m_levels, m_log, m_wlog, m_ppml))
