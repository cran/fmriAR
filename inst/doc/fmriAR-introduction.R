## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#",
  fig.width = 7,
  fig.height = 5
)

## ----simulate-----------------------------------------------------------------
library(fmriAR)
set.seed(42)

n_time <- 240
n_vox <- 60
runs <- rep(1:2, each = n_time / 2)

# Design matrix: intercept + task regressor
X <- cbind(
  intercept = 1,
  task = rep(c(rep(0, 30), rep(1, 30)), length.out = n_time)
)

# Generate voxel time-series with AR(2) noise
phi_true <- matrix(0, nrow = n_vox, ncol = 2)
phi_true[, 1] <- 0.5 + rnorm(n_vox, sd = 0.05)
phi_true[, 2] <- -0.2 + rnorm(n_vox, sd = 0.05)

Y <- matrix(0, n_time, n_vox)
innov <- matrix(rnorm(n_time * n_vox, sd = 1), n_time, n_vox)
for (v in seq_len(n_vox)) {
  for (t in seq_len(n_time)) {
    ar_part <- 0
    if (t > 1) ar_part <- ar_part + phi_true[v, 1] * Y[t - 1, v]
    if (t > 2) ar_part <- ar_part + phi_true[v, 2] * Y[t - 2, v]
    Y[t, v] <- ar_part + innov[t, v]
  }
}

# Add task signal to half the voxels
beta <- c(0, 1.5)
Y_signal <- drop(X %*% beta)
Y[, 1:30] <- sweep(Y[, 1:30], 1, Y_signal, "+")

# Residuals from initial OLS fit
coeff_ols <- qr.solve(X, Y)
resid <- Y - X %*% coeff_ols

## ----fit-ar-------------------------------------------------------------------
plan_ar <- fit_noise(
  resid,
  runs = runs,
  method = "ar",
  p = "auto",
  p_max = 6,
  exact_first = "ar1",
  pooling = "run"
)

plan_ar

## ----whiten-run---------------------------------------------------------------
whitened <- whiten_apply(plan_ar, X, Y, runs = runs)
str(whitened)

## ----gls----------------------------------------------------------------------
Xw <- whitened$X
Yw <- whitened$Y
beta_gls <- qr.solve(Xw, Yw[, 1:5])
beta_gls

## ----whiteness----------------------------------------------------------------
innov_var <- Yw - Xw %*% qr.solve(Xw, Yw)
lag_stats <- apply(innov_var, 2, function(y) {
  ac <- acf(y, plot = FALSE, lag.max = 5)$acf[-1]
  mean(abs(ac))
})
mean(lag_stats)

## ----whiteness-raw------------------------------------------------------------
lag_stats_raw <- apply(resid, 2, function(y) {
  ac <- acf(y, plot = FALSE, lag.max = 5)$acf[-1]
  mean(abs(ac))
})
mean(lag_stats_raw)

## ----whiteness-plot, fig.cap = "Average autocorrelation across voxels before and after whitening"----
avg_acf <- function(mat, lag.max = 20) {
  acf_vals <- sapply(seq_len(ncol(mat)), function(j) {
    stats::acf(mat[, j], plot = FALSE, lag.max = lag.max)$acf
  })
  rowMeans(acf_vals)
}

lag_max <- 20
lags <- seq_len(lag_max)
acf_raw <- avg_acf(resid, lag.max = lag_max)
acf_white <- avg_acf(innov_var, lag.max = lag_max)

ylim <- range(c(acf_raw[-1L], acf_white[-1L], 0))
plot(lags, acf_raw[-1L], type = "h", lwd = 2, col = "#1b9e77",
     ylim = ylim, xlab = "Lag", ylab = "Average autocorrelation",
     main = "Mean autocorrelation across voxels")
lines(lags, acf_white[-1L], type = "h", lwd = 2, col = "#d95f02")
abline(h = 0, col = "grey70", lty = 2)
legend("topright",
       legend = c("Raw residuals", "Whitened innovations"),
       col = c("#1b9e77", "#d95f02"),
       lwd = 2, bty = "n")

## ----multiscale---------------------------------------------------------------
parcels_fine <- rep(1:12, each = n_vox / 12)
parcels_medium <- rep(1:6, each = n_vox / 6)
parcels_coarse <- rep(1:3, each = n_vox / 3)

plan_parcel <- fit_noise(
  resid,
  runs = runs,
  method = "ar",
  p = "auto",
  p_max = 6,
  pooling = "parcel",
  parcels = parcels_fine,
  parcel_sets = list(
    fine = parcels_fine,
    medium = parcels_medium,
    coarse = parcels_coarse
  ),
  multiscale = "pacf_weighted"
)

plan_parcel$order

## ----whiten-parcel------------------------------------------------------------
whitened_parcel <- whiten_apply(plan_parcel, X, Y, runs = runs, parcels = parcels_fine)
length(whitened_parcel$X_by)

## ----arma---------------------------------------------------------------------
plan_arma <- fit_noise(
  resid,
  runs = runs,
  method = "arma",
  p = 2,
  q = 1,
  hr_iter = 1
)

plan_arma$order

## ----arma-whiten--------------------------------------------------------------
whitened_arma <- whiten_apply(plan_arma, X, Y, runs = runs)
whitened_arma$X[1:5, ]

## ----afni-plan, eval = FALSE--------------------------------------------------
# # Example: AFNI-style AR(3) with optional MA(1) term
# roots <- list(a = 0.6, r1 = 0.7, t1 = pi / 6)
# 
# plan_afni <- compat$afni_restricted_plan(
#   resid = resid,
#   runs = runs,
#   parcels = NULL,
#   p = 3L,
#   roots = roots,
#   estimate_ma1 = TRUE
# )
# 
# plan_afni$order
# 
# # Apply whitening just like any other plan
# whitened_afni <- whiten_apply(plan_afni, X, Y, runs = runs)

## ----diagnostics--------------------------------------------------------------
# Autocorrelation diagnostics for whitened residuals
acorr <- acorr_diagnostics(Yw[, 1:3])
acorr

# Sandwich standard errors from whitened residuals
sandwich <- sandwich_from_whitened_resid(whitened$X, whitened$Y[, 1:3])
sandwich$se

