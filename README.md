# Bayesian VAR Models for Asian FX Forecasting

A comparative study of three Bayesian Vector Autoregression (BVAR) models for forecasting Asian currency exchange rates against the US Dollar, with the DXY index as an exogenous driver.

---

## Overview

This project estimates and evaluates three competing BVAR specifications on daily log-returns of INR/USD, CNY/USD, and SGD/USD (plus DXY), comparing their out-of-sample forecast accuracy using both point and probabilistic loss metrics.

| Model | Key Feature |
|---|---|
| **Standard BVAR** | Minnesota prior with hyperparameter optimization via Metropolis-Hastings |
| **BVAR-SV** | Time-varying stochastic volatility (Primiceri-style) for heteroskedastic innovations |
| **BVAR-t** | Student-t distributed innovations for robustness to fat-tailed exchange rate shocks |

---

## Data

- **Series**: INR/USD, CNY/USD, SGD/USD, DXY (US Dollar Index)
- **Source**: Yahoo Finance via `quantmod`
- **Sample**: January 2, 2010 – April 2, 2026
- **Transformation**: Daily log-returns scaled by 100
- **Train/Test split**: 80% training, 20% test (hold-out)

---

## Models

### 1. Standard BVAR (`BVAR` package)
- Minnesota prior with full hyperparameter optimization
- Lag length selected via AIC on training data (`VARselect`)
- 60,000 MCMC draws, 20,000 burn-in, thinning factor of 2
- Metropolis-Hastings with adaptive acceptance rate (target: 20–40%)

### 2. BVAR with Stochastic Volatility (`bvarsv` package)
- Time-varying parameter model with log-variance random walk
- Log-variance propagated forward for multi-step forecasts
- Forecast draws: 2,000 posterior samples with stability caps applied to prevent explosive paths

### 3. BVAR-t (custom Gibbs sampler)
- Multivariate Student-t innovations with fixed degrees of freedom (`ν = 8`)
- Implemented via data augmentation: auxiliary Gamma weights (ω) scale the covariance per observation
- Vectorized Gibbs sampler for efficiency (5,000 iterations, 1,000 burn-in)

---

## Evaluation

### Metrics
- **RMSE** — point forecast accuracy
- **CRPS** (Continuous Ranked Probability Score) — full predictive distribution sharpness and calibration
- **Diebold-Mariano test** — formal test of equal predictive accuracy between models (two-sided, CRPS loss)

### Outputs
- RMSE and CRPS comparison tables across all currencies and models
- DM test statistics and p-values for BVAR vs BVAR-SV and BVAR vs BVAR-t
- Time series plot of CRPS improvement over standard BVAR across the forecast horizon

---

## Repository Structure

```
.
├── main.R              # Full analysis pipeline (data → models → evaluation → plots)
└── README.md
```

---

## Requirements

**R ≥ 4.2** with the following packages:

```r
quantmod, ggplot2, tidyr, dplyr, BVAR, vars, coda,
patchwork, lubridate, scoringRules, forecast, bvarsv, MASS
```

All packages are auto-installed on first run if missing.

---

## Usage

```r
# Clone the repo and run the full pipeline
source("main.R")
```

The script runs sequentially:
1. Downloads and cleans data
2. Produces EDA plots (levels, returns, correlations, densities)
3. Selects lag length via AIC
4. Estimates all three BVAR models
5. Generates forecasts and computes RMSE / CRPS
6. Runs Diebold-Mariano tests
7. Plots model improvement over the forecast horizon

> **Note**: BVAR-SV estimation (`bvar.sv.tvp`) and the BVAR-t Gibbs sampler are computationally intensive. Runtime on a standard laptop is approximately 15–30 minutes depending on hardware.

---

## Key Design Choices

- **Stability caps on BVAR-SV forecasts**: Log-variance is clamped to `[-10, 5]` and individual forecasts to `[-20, 20]` to prevent explosive draws in long-horizon simulation. These are presented as robustness adjustments and should be discussed in any write-up.
- **Uniform lag length**: All three models use the same AIC-selected lag `p` for a fair comparison.
- **CRPS as primary metric**: RMSE only evaluates point forecasts; CRPS rewards well-calibrated uncertainty, which is the more relevant criterion for risk management applications.

---

## Results Interpretation

- Values > 0 in the improvement plot indicate the extension outperforms the standard BVAR on CRPS.
- A statistically significant DM test (p < 0.05) indicates the accuracy difference is not due to chance.
- BVAR-SV is expected to outperform during high-volatility regimes (e.g., COVID-19, 2022 USD surge); BVAR-t is expected to be more robust to isolated large shocks.

---

