# ==============================================================================
# DiD estimation script (R), companion to Winkler et al. (2026):
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
# Requires: fixest. For a version-pinned run, install renv once and keep the renv::use() line below.
# ==============================================================================

# Reproducible environment: pin the exact package version (renv >= 1.0).
# Activates an isolated library; comment out to use your own installed fixest.
renv::use("fixest@0.11.2")

library(fixest)

df <- read.csv("did_lab_panel.csv")
cat(sprintf("rows = %d   sum(y) = %.0f\n", nrow(df), sum(as.numeric(df$y))))

# scale for the levels estimate: treated pre-mean x control post/pre ratio
mean_treat_pre <- mean(df$y[df$treat == 1 & df$post == 0])
mean_ctrl_pre  <- mean(df$y[df$treat == 0 & df$post == 0])
mean_ctrl_post <- mean(df$y[df$treat == 0 & df$post == 1])
scale          <- mean_treat_pre * mean_ctrl_post / mean_ctrl_pre

n_pairs <- length(unique(df$pair))
tcrit   <- qt(0.975, df = n_pairs - 1)

df$log_y <- log1p(df$y)
w_pre    <- tapply(df$y[df$post == 0], df$unit[df$post == 0], mean)
df$w_pre <- as.numeric(w_pre[as.character(df$unit)])

# SEs: CRV1 clustered by matched pair. With 1:1 matching WITHOUT replacement
# (this panel), each unit belongs to exactly one pair, so two-way unit-and-pair
# clustering is identical to pair clustering. If you match WITH replacement
# (controls reused across pairs), units no longer nest in pairs, so cluster
# two-way instead: cluster = ~unit + pair.
m_levels <- feols(y     ~ D | unit + period, data = df, cluster = ~pair)
m_log    <- feols(log_y ~ D | unit + period, data = df, cluster = ~pair)
m_wlog   <- feols(log_y ~ D | unit + period, data = df[df$w_pre > 0, ],
                  weights = ~w_pre, cluster = ~pair)  # zero-weight units drop out anyway
m_ppml   <- fepois(y    ~ D | unit + period, data = df, cluster = ~pair)

# delta-method CI for the scaled levels effect theta = b/scale: the scale is
# estimated from the same sample and strongly negatively correlated with b,
# so a fixed-scale CI would be ~2.5x too wide.
D_dem    <- demean(X = data.frame(D = df$D), f = df[, c("unit", "period")], tol = 1e-8)[, 1]
resid_lv <- resid(m_levels)
denom    <- sum(D_dem^2)
score_pair   <- tapply(D_dem * resid_lv, df$pair, sum)
pm_treat_pre <- tapply(df$y[df$treat == 1 & df$post == 0], df$pair[df$treat == 1 & df$post == 0], mean)
pm_ctrl_pre  <- tapply(df$y[df$treat == 0 & df$post == 0], df$pair[df$treat == 0 & df$post == 0], mean)
pm_ctrl_post <- tapply(df$y[df$treat == 0 & df$post == 1], df$pair[df$treat == 0 & df$post == 1], mean)
theta <- unname(coef(m_levels)["D"]) / scale
psi <- score_pair / (denom * scale) -
  (theta / n_pairs) * ((pm_treat_pre - mean(pm_treat_pre)) / mean(pm_treat_pre) -
                       (pm_ctrl_pre  - mean(pm_ctrl_pre))  / mean(pm_ctrl_pre)  +
                       (pm_ctrl_post - mean(pm_ctrl_post)) / mean(pm_ctrl_post))
# CRV1, fixest-exact: BOTH small-sample factors multiply the VARIANCE (inside the sqrt).
# Nested K = period FE + D (unit FE absorbed within the pair clusters, excluded).
K_dof    <- 1 + length(unique(df$period))
ssc      <- (n_pairs / (n_pairs - 1)) * ((nrow(df) - 1) / (nrow(df) - K_dof))
se_theta <- sqrt(sum(psi^2) * ssc)

ci <- function(est, s) sprintf("[%+.4f%%, %+.4f%%]", 100*(est - tcrit*s), 100*(est + tcrit*s))
cat(sprintf("levels  %+.8f   %s   (delta method)\n", theta, ci(theta, se_theta)))
cat(sprintf("log1y   %+.8f   %s\n", coef(m_log)["D"],  ci(coef(m_log)["D"],  se(m_log)["D"])))
cat(sprintf("wlog    %+.8f   %s\n", coef(m_wlog)["D"], ci(coef(m_wlog)["D"], se(m_wlog)["D"])))
d   <- coef(m_ppml)["D"]
sd_ <- se(m_ppml)["D"]
cat(sprintf("ppml    %+.8f   [%+.4f%%, %+.4f%%]\n", exp(d) - 1,
            100*(exp(d - tcrit*sd_) - 1), 100*(exp(d + tcrit*sd_) - 1)))


# Standard regression table (raw coefficients).
# NOTE: levels is in outcome units here (not the scaled % above), PPML is the
# log coefficient (not exp(b)-1), and SEs use fixest's default convention.
cat("\n")
etable(m_levels, m_log, m_wlog, m_ppml)
