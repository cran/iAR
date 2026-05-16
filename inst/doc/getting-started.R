## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----load-pkg-----------------------------------------------------------------
library(iAR)

## ----gentime------------------------------------------------------------------
set.seed(2847)

# x defaults to NULL — no setup needed
times <- gentime(n = 100)@times

cat("Number of observations:", length(times), "\n")
cat("Time span: [", round(min(times), 1), ",", round(max(times), 1), "]\n")
cat("Mean gap between observations:", round(mean(diff(times)), 2), "\n")

## ----gentime-plot, fig.cap = "Distribution of inter-observation gaps (log scale) for 100 irregular times drawn from an exponential mixture."----
gaps <- diff(times)
hist(log1p(gaps),
     main = "Inter-observation gaps (log1p scale)",
     xlab = "log(1 + gap)", col = "steelblue", border = "white",
     breaks = 20)

## ----iar-sim------------------------------------------------------------------
set.seed(2847)
times <- gentime(n = 100)@times

# Build the model: family = "norm", phi = 0.9
m_iar <- iAR(family = "norm", times = times, coef = 0.9)
m_iar <- sim(m_iar)

head(m_iar@series)

## ----iar-plot, fig.cap = "Simulated iAR series (phi = 0.9). The irregular spacing is visible in the unequal distances along the x-axis."----
plot(m_iar, type = "o",main = "Simulated iAR series (phi = 0.9)",ylab="Values",xlab="Time",pch=20)

## ----iar-prep-----------------------------------------------------------------
y       <- m_iar@series
y_std   <- y / sd(y)           # unit variance
m_iar@series     <- y_std
m_iar@series_esd <- rep(0, length(times))

## ----iar-loglik---------------------------------------------------------------
m_loglik <- loglik(m_iar)
cat("loglik estimate: phi =", round(m_loglik@coef, 4), "\n")
cat("Log-likelihood: ", round(m_loglik@loglik, 4), "\n")

## ----iar-kalman---------------------------------------------------------------
m_kalman <- kalman(m_iar)
cat("kalman estimate: phi =", round(m_kalman@coef, 4), "\n")
cat("Kalman log-lik: ", round(m_kalman@kalmanlik, 4), "\n")

## ----iar-hessian--------------------------------------------------------------
m_h <- iAR(family = "norm", times = times, coef = 0.9, hessian = TRUE)
m_h@series     <- y_std
m_h@series_esd <- rep(0, length(times))
m_h <- loglik(m_h)

# Summary contains the coefficient, SE, p-value, information criteria, and residual diagnostics
summary(m_h)

## ----iar-forecast-------------------------------------------------------------
m_fc <- forecast(m_kalman, tAhead = 10)
cat("Forecast (10 time units ahead):", round(m_fc@forecast, 4), "\n")

## ----iar-plot-forecast, fig.cap = "iAR series with one-step-ahead forecast (blue dot)."----
plot_forecast(m_fc,ylab="Values",xlab="Time",type="o",pch=20)

## ----iar-interp---------------------------------------------------------------
# Introduce a missing value at position 10
m_interp           <- m_kalman
m_interp@series[10] <- NA

m_interp <- interpolation(m_interp)
cat("Interpolated value at position 10:",
    round(m_interp@interpolated_values, 4), "\n")
cat("Original value at position 10    :",
    round(y_std[10], 4), "\n")

## ----ciar-sim-----------------------------------------------------------------
set.seed(2847)
times <- gentime(n = 100)@times

# phi = (0.9, 0): purely real CiAR — same as iAR but more general
m_ciar <- CiAR(times = times, coef = c(0.9, 0))
m_ciar <- sim(m_ciar)

# Standardise
y_c           <- m_ciar@series
y_c_std       <- y_c / sd(y_c)
m_ciar@series     <- y_c_std
m_ciar@series_esd <- rep(0, length(times))

## ----ciar-kalman--------------------------------------------------------------
m_ciar <- kalman(m_ciar)
cat("CiAR estimate: phiR =", round(m_ciar@coef[1], 4),
    " phiI =", round(m_ciar@coef[2], 4), "\n")
cat("Modulus |phi| =", round(sqrt(sum(m_ciar@coef^2)), 4), "\n")

## ----ciar-forecast------------------------------------------------------------
m_ciar <- forecast(m_ciar, tAhead = 10)
cat("CiAR forecast:", round(m_ciar@forecast, 4), "\n")

## ----biar-sim-----------------------------------------------------------------
set.seed(2847)
times <- gentime(n = 100)@times

# coef = c(phiR, phiI), rho = cross-series correlation
m_biar <- BiAR(times = times, coef = c(0.9, 0.3), rho = 0.9)
m_biar <- sim(m_biar)

head(m_biar@series)   # matrix: column 1 = series 1, column 2 = series 2

## ----biar-kalman--------------------------------------------------------------
n       <- length(times)
y_b     <- m_biar@series
y_b_std <- y_b / matrix(apply(y_b, 2, sd), nrow = n, ncol = 2, byrow = TRUE)

m_biar@series     <- y_b_std
m_biar@series_esd <- matrix(0, n, 2)

m_biar <- kalman(m_biar)
cat("BiAR estimate: phiR =", round(m_biar@coef[1], 4),
    " phiI =", round(m_biar@coef[2], 4), "\n")

## ----biar-forecast------------------------------------------------------------
m_biar <- forecast(m_biar, tAhead = 10)
cat("BiAR forecast (series 1):", round(m_biar@forecast[1], 4), "\n")
cat("BiAR forecast (series 2):", round(m_biar@forecast[2], 4), "\n")

## ----agn-load-----------------------------------------------------------------
data(agn)
head(agn)

## ----agn-plot, fig.cap = "AGN MCG-6-30-15 light curve. Gaps in the observation schedule produce the irregular spacing typical of ground-based astronomical surveys."----
plot(agn$t, agn$m, xlab = "Heliocentric JD - 2450000",
     ylab = "Flux [10^-15 erg/s/cm^2/A]",
     main = "AGN MCG-6-30-15 (K band)",type="o",pch=20)

## ----agn-fit------------------------------------------------------------------
# Standardise: centre and scale by SD (measurement errors available)
y_agn   <- agn$m
y_agn_c <- (y_agn - mean(y_agn)) / sd(y_agn)
esd_agn <- agn$merr / sd(y_agn)

m_agn <- iAR(
  family     = "norm",
  times      = agn$t,
  series     = y_agn_c,
  series_esd = esd_agn
)

m_agn <- kalman(m_agn)
cat("Estimated phi:", round(m_agn@coef, 4), "\n")
cat("Kalman log-lik:", round(m_agn@kalmanlik, 4), "\n")

## ----agn-fit-plot, fig.cap = "AGN light curve (standardised) with fitted iAR values overlaid."----
m_agn <- fit(m_agn)
plot_fit(m_agn, xlab = "Heliocentric JD - 2450000",
     ylab = "Flux [10^-15 erg/s/cm^2/A]",type="o",pch=20)

## ----agn-forecast-------------------------------------------------------------
m_agn <- forecast(m_agn, tAhead = 1)
cat("One-step-ahead forecast:", round(m_agn@forecast, 4), "\n")

## ----pair---------------------------------------------------------------------
data(cvnovag)
data(cvnovar)

# Pair the G-band and R-band CV Nova light curves
paired <- pairingits(data1 = cvnovag, data2 = cvnovar, tol = 0.01)

cat("Paired observations:", nrow(paired@paired), "\n")
head(paired@paired)

