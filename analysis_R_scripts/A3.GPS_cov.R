library(mgcv)      # GAM
library(spdep)     # spatial weights (if needed later)
library(tidyverse) # data manipulation & plotting
library(tableone)  # covariate balance tables
library(cobalt)    # love.plot for balance
library(survey)    # for weighted analyses
library(broom)     # tidy model outputs

# =============================================================================
# 1 Load & prepare data -----------
# =============================================================================
base_dir <- "C:/Users/anapa/OneDrive - University of Cambridge/PhD - Research/P2 - Spatial Analysis/R Analysis/3rdAnalysis"
setwd(base_dir)

hex_data <- read_rds("Data/Lag/Hex_lag13.rds")
dir.create("Results/GPS", showWarnings = FALSE)
# =============================================================================
# 2. Define variable groups ----
# =============================================================================
id_var <- "H_ID"

treatment <- c("N1_8Area_log_z")

confounders <- c("H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z", 
                "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z", "H_CrowdP_z")

prognostics <- c("H_GreenPN_z", "H_PublicP_z", "H_ActiveP_log_z")

pca_lags <- c("lag13_covPC1", "lag13_covPC2", "lag13_covPC3")

# =============================================================================
# 3. Prepare data ----
# =============================================================================
df <- hex_data %>%
  mutate(
    # create quintiles of treatment for overlap diagnostics
    treat_q = ntile(N1_8Area_log_z, 5))

# =============================================================================
# 4.   Build GPS formulas ----
# =============================================================================
all_gps_vars <- c(confounders, prognostics, pca_lags)

gps_f_lin    <- as.formula(
  paste("N1_8Area_log_z ~", paste(all_gps_vars, collapse = " + "))
)
gps_f_smooth <- as.formula(
  paste("N1_8Area_log_z ~", paste(sprintf("s(%s)", all_gps_vars), collapse = " + ")
  ))

# =============================================================================
# 4.   Estimate GPS Models (lm and gam) ----
# =============================================================================
# 4a. Linear GPS
lm_gps <- lm(gps_f_lin, data = df)
df <- df %>%
  mutate(
    mu_lm  = predict(lm_gps, newdata = df),
    sd_lm  = sigma(lm_gps),
    gps_lm = dnorm(N1_8Area_log_z, mean = mu_lm, sd = sd_lm)
  )

# 4b. GAM GPS
gam_gps <- gam(gps_f_smooth, data = df)
resid_sd <- sd(residuals(gam_gps))
df <- df %>%
  mutate(
    mu_gam  = predict(gam_gps, newdata = df),
    sd_gam  = resid_sd,
    gps_gam = dnorm(N1_8Area_log_z, mean = mu_gam, sd = sd_gam)
  )

# =============================================================================
# 5.   Compare fits ----
# =============================================================================

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
metrics <- tibble(
  model   = c("LM", "GAM"),
  AIC     = AIC(lm_gps, gam_gps)$AIC,
  BIC     = BIC(lm_gps, gam_gps)$BIC,
  RMSE    = c(rmse(df$N1_8Area_log_z, df$mu_lm), rmse(df$N1_8Area_log_z, df$mu_gam)),
  LogLik  = c(logLik(lm_gps)[1], logLik(gam_gps)[1]),
  CorrObs = c(cor(df$N1_8Area_log_z, df$mu_lm), cor(df$N1_8Area_log_z, df$mu_gam))
)
print(metrics)
write.csv(metrics,  "Results/GPS/gpscov_comparefit_lm_gam.csv", row.names = FALSE)

# =============================================================================
# 6.   Compute stabilized GPS weights ----
# =============================================================================

df <- df %>%
  mutate(
    f_t     = approx(density(N1_8Area_log_z)$x, density(N1_8Area_log_z)$y, xout = N1_8Area_log_z)$y,
    sw_raw  = f_t / gps_gam
  )

# =============================================================================
# 7.   Trimming & winsorizing ----
# =============================================================================

# Choose cutoff (e.g. 99th percentile of sw_raw)
cut99 <- quantile(df$sw_raw, .99, na.rm=TRUE)

# A) Trimmed dataset + normalized weights
df_trim   <- df %>%
  filter(sw_raw <= cut99) %>%
  mutate(sw_trim = sw_raw / mean(sw_raw, na.rm=TRUE))

# B) Winsorized weights on full data + normalization
df_winsor <- df %>%
  mutate(
    sw_win = pmin(sw_raw, cut99),
    sw_win = sw_win / mean(sw_win, na.rm=TRUE)
  )

# =============================================================================
# 8.   Covariate balance checks ----
# =============================================================================

run_balance <- function(data, weight_col = NULL) {
  message("\n>>> Balance ", if (is.null(weight_col)) "unweighted" else weight_col, " <<<")
  vars <- c(confounders, prognostics, pca_lags)
  frm  <- as.formula(paste("treatment ~", paste(vars, collapse = " + ")))
  
  if (is.null(weight_col)) {
    bal <- bal.tab(x = frm, data = data, m.threshold = 0.1, quick = FALSE)
  } else {
    bal <- bal.tab(x = frm, data = data,
                   weights  = data[[weight_col]],
                   estimand = "ATE",
                   method   = "weighting",
                   m.threshold = 0.1,
                   quick    = FALSE)
    love.plot(bal, stat = "mean.diffs", threshold = 0.1)
  }
  print(bal)
}

# 8a) Unweighted
run_balance(df)
# 8b) Raw GPS weights
run_balance(df,       weight_col = "sw_raw")
# 8c) Trimmed stabilized weights
run_balance(df_trim,  weight_col = "sw_trim")
# 8d) Winsorized stabilized weights
run_balance(df_winsor,weight_col = "sw_win")

# =============================================================================
# 9.   Overlap Plots ----
# =============================================================================
# Raw
gg_raw <- df %>%
  ggplot(aes(x = gps_gam,
             weight = sw_raw,
             fill   = factor(treat_q))) +
  geom_density(alpha = .5) +
  labs(title = "Weighted GAM GPS density by treatment quintile (raw)",
       x     = "GPS (GAM)",
       fill  = "Quintile") +
  theme_minimal()
ggsave("Results/GPS/overlap_rawcov.png", gg_raw, width = 6, height = 4)

# Winsorized
gg_winsor <- df_winsor %>%
  ggplot(aes(x      = gps_gam,
             weight = sw_win,
             fill   = factor(treat_q))) +
  geom_density(alpha = .5) +
  labs(title = "Weighted GAM GPS density by treatment quintile (winsorized)",
       x = "GPS (GAM)", fill = "Quintile") +
  theme_minimal()
# Save plot
ggsave("Results/GPS/overlap_winsorcov.png", gg_winsor, width = 6, height = 4)

# Trimmed
gg_trim <- df_trim %>%
  ggplot(aes(x      = gps_gam,
             weight = sw_trim,
             fill   = factor(treat_q))) +
  geom_density(alpha = .5) +
  labs(title = "Weighted GAM GPS density by treatment quintile (trimmed)",
       x = "GPS (GAM)", fill = "Quintile") +
  theme_minimal()
# Save plot
ggsave("Results/GPS/overlap_trimcov.png", gg_trim, width = 6, height = 4)

# =============================================================================
# 10.   Weight distribution diagnostics ----
# =============================================================================

# Histogram of raw stabilized weights
gg_hist_raw <- df %>%
  ggplot(aes(x = sw_raw)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of raw stabilized GPS weights", x = "sw_raw")
# Save plot
ggsave("Results/GPS/hist_sw_rawcov.png", gg_hist_raw, width = 6, height = 4)

# Histogram of Trimmed stabilized weights
gg_hist_trim <- df_trim %>%
  ggplot(aes(x = sw_trim)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of trimmed stabilized GPS weights", x = "sw_trim")
# Save plot
ggsave("Results/GPS/hist_sw_trimcov.png", gg_hist_trim, width = 6, height = 4)

# Histogram of Winsorized stabilized weights
gg_hist_win <- df_winsor %>%
  ggplot(aes(x = sw_win)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of winsorized stabilized GPS weights", x = "sw_win")
# Save plot
ggsave("Results/GPS/hist_sw_wincov.png", gg_hist_win, width = 6, height = 4)

quant_raw <- df %>% summarise(across(sw_raw, ~quantile(.x, probs = c(0.01,0.05,0.5,0.95,0.99)))) %>%
  pivot_longer(everything(), names_to = "weight", values_to = "quantiles") %>%
  unnest_wider(quantiles, names_sep = "_")
write_csv(quant_raw, "Results/GPS/quantiles_sw_rawcov.csv")

quant_trim <- df_trim %>% summarise(across(sw_trim, ~quantile(.x, probs = c(0.01,0.05,0.5,0.95,0.99)))) %>%
  pivot_longer(everything(), names_to = "weight", values_to = "quantiles") %>%
  unnest_wider(quantiles, names_sep = "_")
write_csv(quant_trim, "Results/GPS/quantiles_sw_trimcov.csv")

quant_win <- df_winsor %>% summarise(across(sw_win, ~quantile(.x, probs = c(0.01,0.05,0.5,0.95,0.99)))) %>%
  pivot_longer(everything(), names_to = "weight", values_to = "quantiles") %>%
  unnest_wider(quantiles, names_sep = "_")
write_csv(quant_win, "Results/GPS/quantiles_sw_wincov.csv")

# =============================================================================
# 11.   Standardize GPS ----
# =============================================================================

# Compute skewness of GPS 
gps_skew   <- psych::skew(df$gps_gam, na.rm = TRUE)
skew_thresh <- 1

# If negative skew, precompute max for reflection
gps_max <- max(df$gps_gam, na.rm = TRUE)

# 3. Branch on the scalar skewness, then mutate
if (gps_skew >  skew_thresh) {
  # strong right skew → simple log
  message("Right‐skew → using log(gps + 1)")
  df <- df %>%
    mutate(gps_log = log(gps_gam + 1))
  
} else if (gps_skew < -skew_thresh) {
  # strong left skew → reflect, log, then negate
  message("Left‐skew → using -log(max(gps)+1 - gps)")
  df <- df %>%
    mutate(
      gps_rlog = -log(gps_max + 1 - gps_gam))
} 
"gps_log" %in% names(df)
"gps_rlog" %in% names(df)

# 4. Now standardize the new gps_log into gps_log_z
df <- df %>%
  mutate(
    gps_gam_z = as.numeric(scale(gps_gam))
  )
df_trim <- df_trim %>%
  mutate(
    gps_gam_z = as.numeric(scale(gps_gam))
  )
df_winsor <- df_winsor %>%
  mutate(
    gps_gam_z = as.numeric(scale(gps_gam))
  )

# =============================================================================
# 12.   Save all results for later use & reporting ----
# =============================================================================

# 9a) Data
write_rds(df,       "Results/GPS/gps_weightscov.rds")
write_rds(df_trim,  "Results/GPS/gps_trimmedcov.rds")
write_rds(df_winsor,"Results/GPS/gps_winsorizedcov.rds")

# 9b) Models & metrics
write_rds(lm_gps,  "Results/GPS/lm_gpscov.rds")
write_rds(gam_gps, "Results/GPS/gam_gpscov.rds")

# 9c) Balance tables
df_bal_unw <- as.data.frame(bal.tab(treatment ~ ., data = df, weights = rep(1,nrow(df)))$Balance)
write_csv(df_bal_unw, "Results/GPS/balance_unweightedcov.csv")
df_bal_trim <- as.data.frame(bal.tab(treatment ~ ., data = df_trim, weights = df_trim$sw_trim)$Balance)
write_csv(df_bal_trim, "Results/GPS/balance_trimmedcov.csv")
df_bal_win  <- as.data.frame(bal.tab(treatment ~ ., data = df_winsor, weights = df_winsor$sw_win)$Balance)
write_csv(df_bal_win,  "Results/GPS/balance_winsorizedcov.csv")

