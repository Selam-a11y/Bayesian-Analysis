# =============================================================================
# 0. SETUP
# =============================================================================

required_packages <- c(
  "quantmod", "ggplot2", "tidyr", "dplyr", "BVAR",
  "vars", "coda", "patchwork", "lubridate",
  "scoringRules", "forecast", "bvarsv", "MASS"
)

installed <- installed.packages()[, "Package"]
new_pkg <- required_packages[!(required_packages %in% installed)]
if (length(new_pkg) > 0) install.packages(new_pkg)

invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)


# =============================================================================
# 1. DATA DOWNLOAD AND CLEANING
# =============================================================================

symbols <- c("INR=X", "CNY=X", "SGD=X", "DX-Y.NYB")

getSymbols(
  symbols,
  src = "yahoo",
  from = "2010-01-02",
  to   = "2026-04-02",
  auto.assign = TRUE
)

data_fx <- data.frame(
  Date = index(`INR=X`),
  INR  = as.numeric(`INR=X`[, 4]),
  CNY  = as.numeric(`CNY=X`[, 4]),
  SGD  = as.numeric(`SGD=X`[, 4])
)

data_dxy <- data.frame(
  Date = index(`DX-Y.NYB`),
  DXY  = as.numeric(`DX-Y.NYB`[, 4])
)

data_full <- inner_join(data_fx, data_dxy, by = "Date") |> 
  na.omit()

data_ret <- data_full |>
  mutate(across(c(INR, CNY, SGD, DXY), ~ 100 * log(.x / lag(.x)))) |>
  na.omit()

model_data <- as.matrix(data_ret[, c("INR", "CNY", "SGD", "DXY")])
dates_ret <- data_ret$Date


# =============================================================================
# 2. EXPLORATORY DATA ANALYSIS
# =============================================================================

p_levels <- data_full |>
  pivot_longer(-Date, names_to = "Series", values_to = "Value") |>
  ggplot(aes(Date, Value)) +
  geom_line() +
  facet_wrap(~Series, scales = "free_y") +
  theme_minimal() +
  labs(title = "Exchange Rate Levels")

p_returns <- data_ret |>
  pivot_longer(-Date, names_to = "Series", values_to = "Return") |>
  ggplot(aes(Date, Return)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~Series, scales = "free_y") +
  theme_minimal() +
  labs(title = "Daily Log Returns, Scaled by 100")

cor_matrix <- cor(model_data)

p_corr <- as.data.frame(as.table(cor_matrix)) |>
  ggplot(aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2))) +
  scale_fill_gradient2(
    low = "darkred",
    mid = "white",
    high = "darkblue",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(title = "Return Correlation Matrix", x = "", y = "")

p_density <- as.data.frame(model_data) |>
  pivot_longer(everything(), names_to = "Series", values_to = "Return") |>
  ggplot(aes(Return, fill = Series)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Return Distributions")

print(p_levels)
print(p_returns)
print(p_corr)
print(p_density)


# =============================================================================
# 3. TRAIN-TEST SPLIT AND LAG SELECTION
# =============================================================================

n_total <- nrow(model_data)
n_train <- floor(0.8 * n_total)
h <- n_total - n_train

train_data <- model_data[1:n_train, ]
test_data  <- model_data[(n_train + 1):n_total, ]

lag_select <- VARselect(train_data, lag.max = 5, type = "const")
optimal_lag <- as.integer(lag_select$selection["AIC(n)"])

# For fair comparison, use the same lag everywhere.
p <- optimal_lag

cat("Selected lag length:", p, "\n")



# =============================================================================
# 5. STANDARD BVAR
# =============================================================================

priors_bvar <- bv_priors(
  hyper = "full",
  mn = bv_minnesota(
    lambda = bv_lambda(mode = 0.2, sd = 0.4),
    alpha  = bv_alpha(mode = 2, sd = 0.15),
    #psi    = bv_psi(mode = 0.5, sd = 0.2)
  )
)

mh_setup <- bv_metropolis(
  scale_hess = 0.001,
  adjust_acc = TRUE,
  acc_lower = 0.20,
  acc_upper = 0.40,
  adjust_burn = 0.75
)

fit_bvar <- bvar(
  data   = train_data,
  lags   = p,
  n_draw = 60000,
  n_burn = 20000,
  thin   = 2,
  priors = priors_bvar,
  mh     = mh_setup
)

summary(fit_bvar)
shocks<-irf(fit_bvar, horizon=10)
plot(shocks)

pred_bvar <- predict(fit_bvar, horizon = h)


fc_bvar_draws <- pred_bvar$fcast
fc_bvar_draws <- aperm(fc_bvar_draws, c(2, 3, 1))

fc_bvar_mean <- apply(fc_bvar_draws, c(1, 2), mean)
colnames(fc_bvar_mean) <- colnames(train_data)

# =============================================================================
# BVAR. DIAGNOSTICS

# =============================================================================

bvar_resids <- residuals(fit_bvar)

par(mfrow = c(2, 2))
for (j in seq_len(ncol(bvar_resids))) {
  acf(bvar_resids[, j], main = paste("Residual ACF:", colnames(train_data)[j]))
}
par(mfrow = c(1, 1))

plot(fit_bvar)

# =============================================================================
# 6. BVAR-SV
# =============================================================================

fit_sv <- bvar.sv.tvp(
  train_data,
  p       = p,
  tau     = floor(0.1 * nrow(train_data)),
  nrep    = 5000,
  nburn   = 1000,
  thinfac = 5
)
message("=== BVAR-SV estimation complete ===")

# Important:
# In final project, explain clearly how forecasts are generated.
# Avoid arbitrary caps unless you present them as a robustness adjustment. !!!!!!

# =============================================================================
#  BVAR-SV DIAGNOSTICS
# =============================================================================


# Check the volatility draws for the first variable
# fit_sv$H.draws contains the stochastic volatility draws
plot.ts(fit_sv$H.draws[1, 1, ], main="Trace Plot: Volatility of Var 1", ylab="Value")



T_sv     <- dim(fit_sv$H.draws)[3]
n_stored <- dim(fit_sv$Beta.draws)[1]
n_sims   <- min(2000L, n_stored)
draw_idx <- sample(n_stored, n_sims, replace = FALSE)
k <- ncol(train_data)
fc_sv_draws <- array(NA_real_, dim = c(h, k, n_sims))

message("=== Pre-computing SV innovation SDs from H_draws ===")
sv_innov_sd_all <- matrix(NA_real_, nrow = n_stored, ncol = k)

for (d in seq_len(n_stored)) {
  for (j in seq_len(k)) {
    h_path <- fit_sv$H.draws[j, , d] # variable j, all time, draw d
    # full time path of log-variance
    diffs  <- diff(h_path)                     # first differences ~ N(0, sigma_eta^2)
    sv_innov_sd_all[d, j] <- max(sd(diffs), 1e-8)   # guard against zero
  }
}





message("=== SV innovation SDs computed ===")
message(sprintf("Mean SV innov sd across draws and variables: %.6f",
                mean(sv_innov_sd_all)))





# =============================================================================
#  BVAR-SV FORECAST LOOP 
# =============================================================================

# 1. FIX THE SD CALCULATION (Pre-loop)
# Use the median SD across draws to avoid being spiked by outliers
global_sv_sd <- apply(sv_innov_sd_all, 2, median) 

for (s in seq_len(n_sims)) {
  d <- draw_idx[s]
  
  # Coefficients
  beta_vec <- fit_sv$Beta.draws[d, , T_sv]
  B_draw   <- matrix(beta_vec, nrow = k, ncol = k * optimal_lag + 1, byrow = FALSE)
  
  # Log-variance initialization
  log_var_state <- fit_sv$H.draws[, T_sv, d] 
  
  # Lag initialization
  lag_state <- train_data[(nrow(train_data) - optimal_lag + 1):nrow(train_data), ]
  lag_state <- matrix(lag_state, nrow = optimal_lag, ncol = k)
  
  for (tt in seq_len(h)) {
    
    # A. PROPAGATE WITH SAFETY
    # Use the pre-computed SD but cap the random walk movement
    innov <- rnorm(k, mean = 0, sd = global_sv_sd)
    log_var_state <- log_var_state + innov
    
    # CRITICAL: Tighten the cap. 
    # If your data is daily % returns, log_var shouldn't really exceed 5 (exp(0.5*5) = 12% daily SD)
    log_var_state <- pmax(pmin(log_var_state, 5), -10) 
    sigma_draw    <- exp(0.5 * log_var_state)
    
    # B. GENERATE PREDICTION
    x_vec <- c(1, as.numeric(t(lag_state[optimal_lag:1, , drop = FALSE])))
    y_hat <- as.numeric(B_draw %*% x_vec)
    
    # C. ADD NOISE & VALIDATE
    noise <- rnorm(k, mean = 0, sd = sigma_draw)
    y_new <- y_hat + noise
    
    # EMERGENCY CAP: If y_new is physically impossible (e.g. > 50% move in one day)
    # This stops the feedback loop from exploding the next lag_state
    y_new <- pmax(pmin(y_new, 20), -20) 
    
    fc_sv_draws[tt, , s] <- y_new
    
    # D. ROLL LAGS
    lag_state <- rbind(lag_state[-1, , drop = FALSE], matrix(y_new, nrow = 1, ncol = k))
  }
}


# =============================================================================
# 7. BVAR-t MODEL
# =============================================================================

# This section is acceptable as an extension, but should be clearly separated.
# Present it as a robustness model for fat-tailed exchange-rate innovations.

run_bvar_t <- function(Y, p = 2, nu = 8, n_iter = 5000, burn = 1000) {
  
  # 1. Setup
  T <- nrow(Y)
  k <- ncol(Y)
  
  # Create lags
  embed_Y <- embed(Y, p + 1)
  Y_dep <- embed_Y[, 1:k]
  X <- embed_Y[, -(1:k)]
  X <- cbind(1, X) # Intercept
  
  T_eff <- nrow(Y_dep)
  K <- ncol(X)
  
  # 2. Priors
  B0 <- matrix(0, K, k)
  V0 <- diag(10, K)
  inv_V0 <- solve(V0)
  S0 <- diag(k)
  nu0 <- k + 2
  
  # Initial values
  B <- matrix(0, K, k)
  Sigma <- diag(k)
  omega <- rep(1, T_eff)
  
  # Storage
  B_draws <- array(NA, c(n_iter - burn, K, k))
  Sigma_draws <- array(NA, c(n_iter - burn, k, k))
  
  # 3. Optimized Gibbs Sampler
  for (iter in 1:n_iter) {
    
    # SPEED FIX 1: Avoid diag(omega). Multiply X by omega directly.
    # This scales each observation's influence by its "t-weight"
    XW <- X * omega 
    
    # Update B
    V_post <- solve(t(XW) %*% X + inv_V0)
    B_post <- V_post %*% (t(XW) %*% Y_dep + inv_V0 %*% B0)
    B <- B_post + t(chol(V_post)) %*% matrix(rnorm(K * k), K, k)
    
    # Update Sigma
    E <- Y_dep - X %*% B
    # Use element-wise multiplication for S_post to keep speed up
    S_post <- S0 + t(E * omega) %*% E 
    Sigma <- solve(rWishart(1, T_eff + nu0, solve(S_post))[,,1])
    Sigma_inv <- solve(Sigma)
    
    # SPEED FIX 2: Vectorized Omega sampling
    # Replaces the internal 'for (tt in 1:T_eff)' loop
    quad_all <- rowSums((E %*% Sigma_inv) * E)
    shape_all <- (nu + k) / 2
    rate_all  <- (nu + quad_all) / 2
    
    omega <- rgamma(T_eff, shape = shape_all, rate = rate_all)
    
    # Store results
    if (iter > burn) {
      idx <- iter - burn
      B_draws[idx, , ] <- B
      Sigma_draws[idx, , ] <- Sigma
    }
    
    # Progress indicator
    if (iter %% 500 == 0) message("Iteration: ", iter)
  }
  
  return(list(
    B_draws = B_draws,
    Sigma_draws = Sigma_draws,
    p = p,
    nu = nu
  ))
}

# Run the model
fit_bvart <- run_bvar_t(train_data, p = p)

#diagnostics of bvart

library(coda)
# Check the first coefficient (Intercept for variable 1)
b_chain <- mcmc(fit_bvart$B_draws[, 1, 1])
plot(b_chain)
autocorr.plot(b_chain)
effectiveSize(b_chain)



#forecast for bvar-t


forecast_bvar_t <- function(fit, Y, h = 12) {
  # 1. Setup
  B_draws     <- fit$B_draws
  Sigma_draws <- fit$Sigma_draws
  p           <- fit$p
  nu          <- fit$nu
  n_draws     <- dim(B_draws)[1]
  k           <- ncol(Y)
  
  # Storage for all simulated paths [draws, horizon, variables]
  forecast_storage <- array(NA, c(n_draws, h, k))
  
  # 2. Iterate over each posterior draw
  for (i in 1:n_draws) {
    B_i     <- B_draws[i, , ]     # (K x k) matrix
    Sigma_i <- Sigma_draws[i, , ] # (k x k) matrix
    
    # Get the last 'p' observations to start the forecast
    current_Y_lagged <- Y[(nrow(Y) - p + 1):nrow(Y), ]
    
    # Forecast h steps ahead for this draw
    temp_forecast <- matrix(NA, h, k)
    
    for (step in 1:h) {
      # Prepare predictors (Intercept + Lags)
      # We reverse the lags because embed() and your X matrix use [t-1, t-2...]
      lags_vec <- as.vector(t(current_Y_lagged[p:1, ]))
      x_t      <- c(1, lags_vec)
      
      # Step 1: Calculate the mean prediction
      y_mean <- x_t %*% B_i
      
      # Step 2: Generate t-distributed error
      # Scale by Sigma and the heavy-tail weight (from Inverse Gamma)
      w_t   <- rgamma(1, shape = nu/2, rate = nu/2)
      error <- t(chol(Sigma_i / w_t)) %*% rnorm(k)
      
      # Step 3: Combine and store
      y_next <- y_mean + t(error)
      temp_forecast[step, ] <- y_next
      
      # Update lagged data: drop oldest, add newest
      current_Y_lagged <- rbind(current_Y_lagged[-1, ], y_next)
    }
    
    forecast_storage[i, , ] <- temp_forecast
  }
  
  # 3. Summarize (Mean and Quantiles)
  forecast_mean <- apply(forecast_storage, c(2, 3), mean)
  forecast_q05  <- apply(forecast_storage, c(2, 3), quantile, probs = 0.05)
  forecast_q95  <- apply(forecast_storage, c(2, 3), quantile, probs = 0.95)
  
  return(list(
    mean = forecast_mean,
    lower = forecast_q05,
    upper = forecast_q95,
    all_draws = forecast_storage
  ))
}



# Generate forecasts from the BVAR-t model
fc_t_output <- forecast_bvar_t(fit_bvart, train_data, h = h)

# Extract the Mean (for RMSE) and the Draws (for CRPS)
fc_t_mean <- fc_t_output$mean
# Reorganizing draws to match your SV structure: [horizon, variable, draws]
fc_t_draws <- aperm(fc_t_output$all_draws, c(2, 3, 1))



# =============================================================================
# . FORECAST COMPARISON RMSE AND CRPS BVAR & BVARSV
# =============================================================================

#test_data <- model_data[(n_train + 1):n_total, ]

fc_sv_mean <- apply(fc_sv_draws, c(1, 2), mean)
colnames(fc_sv_mean) <- colnames(train_data)

rmse_fn <- function(pred, act) sqrt(mean((pred - act)^2))
rmse_results <- data.frame(
  Currency = colnames(train_data),
  BVAR     = sapply(seq_len(k), function(j) rmse_fn(fc_bvar_mean[, j], test_data[, j])),
  BVAR_SV  = sapply(seq_len(k), function(j) rmse_fn(fc_sv_mean[,   j], test_data[, j]))
)

message("=== RMSE Results ===")
print(rmse_results)



# Function to calculate CRPS for a specific currency and horizon
calc_crps <- function(draws, actual) {
  # draws: a vector of posterior predictive draws for one point in time
  # actual: the single actual observed value
  return(crps_sample(actual, draws))
}

# Example for the first currency at the first forecast step (h=1)
# BVAR_SV CRPS
sv_crps_val <- calc_crps(fc_sv_draws[1, 1, ], test_data[1, 1])

# You would do the same for your BVAR draws to compare them.

library(scoringRules)

# Initialize storage for CRPS scores
crps_sv <- matrix(NA, nrow = h, ncol = k)
crps_bvar <- matrix(NA, nrow = h, ncol = k)

for (j in seq_len(k)) {
  for (t in seq_len(h)) {
    # CRPS for BVAR-SV
    crps_sv[t, j] <- crps_sample(test_data[t, j], fc_sv_draws[t, j, ])

    # CRPS for BVAR (Assuming you have draws saved similarly)
    crps_bvar[t, j] <- crps_sample(test_data[t, j], fc_bvar_draws[t, j, ])
  }
}

# Aggregate results (Mean CRPS across the forecast horizon)
crps_comparison <- data.frame(
  Currency = colnames(train_data),
  CRPS_BVAR = colMeans(crps_bvar),
  CRPS_BVAR_SV = colMeans(crps_sv)
)

print(crps_comparison)


# =============================================================================
# . FORECAST COMPARISON RMSE AND CRPS BVAR & BVARSV & bvart
# =============================================================================
rmse_results <- data.frame(
  Currency = colnames(train_data),
  BVAR     = sapply(seq_len(k), function(j) rmse_fn(fc_bvar_mean[, j], test_data[, j])),
  BVAR_SV  = sapply(seq_len(k), function(j) rmse_fn(fc_sv_mean[, j], test_data[, j])),
  BVAR_T   = sapply(seq_len(k), function(j) rmse_fn(fc_t_mean[, j], test_data[, j]))
)

message("=== Updated RMSE Results (including BVAR-t) ===")
print(rmse_results)


# Initialize storage for BVAR-t CRPS
crps_bvart <- matrix(NA, nrow = h, ncol = k)

for (j in seq_len(k)) {
  for (t in seq_len(h)) {
    # CRPS for BVAR-SV
    crps_sv[t, j] <- crps_sample(test_data[t, j], fc_sv_draws[t, j, ])

    # CRPS for standard BVAR
    crps_bvar[t, j] <- crps_sample(test_data[t, j], fc_bvar_draws[t, j, ])

    # CRPS for BVAR-t
    crps_bvart[t, j] <- crps_sample(test_data[t, j], fc_t_draws[t, j, ])
  }
}

# Aggregate updated results
crps_comparison <- data.frame(
  Currency     = colnames(train_data),
  CRPS_BVAR    = colMeans(crps_bvar),
  CRPS_BVAR_SV = colMeans(crps_sv),
  CRPS_BVAR_T  = colMeans(crps_bvart)
)

message("=== Updated CRPS Results (including BVAR-t) ===")
print(crps_comparison)

# =============================================================================
# 9. COMBINED 
# =============================================================================

# 1. Prepare RMSE data
rmse_vec_bvar   <- sapply(seq_len(k), function(j) rmse_fn(fc_bvar_mean[, j], test_data[, j]))
rmse_vec_sv     <- sapply(seq_len(k), function(j) rmse_fn(fc_sv_mean[, j], test_data[, j]))
rmse_vec_bvart  <- sapply(seq_len(k), function(j) rmse_fn(fc_t_mean[, j], test_data[, j]))

# 2. Prepare CRPS data (Mean across horizon h)
crps_vec_bvar   <- colMeans(crps_bvar)
crps_vec_sv     <- colMeans(crps_sv)
crps_vec_bvart  <- colMeans(crps_bvart)

# 3. Compile the Final Comparison Table
comparison_table <- data.frame(
  Currency = colnames(train_data),
  
  # RMSE Columns (Lower is better)
  RMSE_BVAR   = round(rmse_vec_bvar, 4),
  RMSE_SV     = round(rmse_vec_sv, 4),
  RMSE_BVAR_T = round(rmse_vec_bvart, 4),
  
  # CRPS Columns (Lower is better)
  CRPS_BVAR   = round(crps_vec_bvar, 4),
  CRPS_SV     = round(crps_vec_sv, 4),
  CRPS_BVAR_T = round(crps_vec_bvart, 4)
)

# 4. Display results
print("--- Final Model Comparison Table ---")
print(comparison_table)


# =============================================================================
# 10. DIEBOLD-MARIANO TEST (FOR CRPS LOSS)
# =============================================================================

library(forecast)

# Define a function to run DM tests across all currencies
run_dm_comparison <- function(loss_matrix_1, loss_matrix_2, model_names = c("M1", "M2")) {
  
  currencies <- colnames(train_data)
  results_list <- list()
  
  for (j in seq_along(currencies)) {
    # The DM test requires the time series of losses
    # We use the CRPS values calculated at each time step 't'
    loss_diff <- loss_matrix_1[, j] - loss_matrix_2[, j]
    
    # dm.test(e1, e2, alternative, h, power)
    # e1 and e2 here are the loss series (CRPS)
    dm_result <- dm.test(loss_matrix_1[, j], 
                         loss_matrix_2[, j], 
                         alternative = "two.sided", 
                         h = 1, # Forecast horizon for the test
                         power = 1)
    
    results_list[[currencies[j]]] <- data.frame(
      Currency = currencies[j],
      Comparison = paste(model_names[1], "vs", model_names[2]),
      DM_Stat = round(dm_result$statistic, 4),
      P_Value = round(dm_result$p.value, 4),
      Significant = ifelse(dm_result$p.value < 0.05, "Yes*", "No")
    )
  }
  
  return(do.call(rbind, results_list))
}

# 1. Compare Standard BVAR vs BVAR-SV
dm_sv <- run_dm_comparison(crps_bvar, crps_sv, c("BVAR", "BVAR-SV"))

# 2. Compare Standard BVAR vs BVAR-t
dm_t <- run_dm_comparison(crps_bvar, crps_bvart, c("BVAR", "BVAR-T"))

# Display Results
cat("\n--- Diebold-Mariano Test: BVAR vs BVAR-SV ---\n")
print(dm_sv)

cat("\n--- Diebold-Mariano Test: BVAR vs BVAR-T ---\n")
print(dm_t)


# =============================================================================
# VISUAL COMPARISON: SV vs T-DISTRIBUTION
# =============================================================================

# 1. Ensure the data frame has both differentials
# We calculate (BVAR_Loss - Model_Loss). 
# Positive values = The extension is improving upon the standard BVAR.
plot_comp_series <- data.frame(
  Date    = dates_ret[(n_train + 1):n_total],
  SV_Gain = rowMeans(crps_bvar - crps_sv),
  T_Gain  = rowMeans(crps_bvar - crps_bvart)
)

# 2. Pivot to long format for easier ggplot mapping
plot_long <- plot_comp_series %>%
  pivot_longer(cols = c(SV_Gain, T_Gain), 
               names_to = "Model", 
               values_to = "Improvement")

# 3. Create the Comparison Plot
ggplot(plot_long, aes(x = Date, y = Improvement, color = Model)) +
  geom_line(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Adding a smoothed trend line to see the 'average' regime shift
  geom_smooth(method = "gam", se = FALSE, size = 1) + 
  scale_color_manual(values = c("SV_Gain" = "#2c7fb8", "T_Gain" = "#31a354"),
                     labels = c("BVAR-SV vs BVAR", "BVAR-t vs BVAR")) +
  labs(title = "Model Improvement over Standard BVAR (CRPS)",
       subtitle = "Values > 0 indicate the extension is outperforming the benchmark",
       y = "Reduction in CRPS ",
       x = "Forecast Horizon Date",
       color = "Model Variant") +
  theme_minimal() +
  theme(legend.position = "bottom")

  

