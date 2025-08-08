###############################################################################
## Subset LM tests for spatial‐lag vs. spatial‐error ##########################
###############################################################################
# 0. Load libraries ----------------------------------------------------------
library(sf)
library(spdep)
library(spatialreg)
library(dplyr)

# 1. Read in your data & weights --------------------------------------------
# Set working directory
base_dir <- "C:/Users/anapa/OneDrive - University of Cambridge/PhD - Research/P2 - Spatial Analysis/R Analysis/3rdAnalysis"
setwd(base_dir)
# Load data
Hex_data       <- read_rds("Data/Mod/Hex_data.rds")
list_IDW8400_13 <- readRDS("Results/Lags/idw8400_13w_filt.rds")
dir.create("Results/Lagrange", showWarnings = FALSE)
colnames(Hex_data)
# 2. Prepare the weighted dataframe (same as before) ------------------------
oc          <- "avgAF_z"
treat       <- "N1_8Area_log_z"
confounders <- c("H_Density_log_z","H_Min65P_z","H_WhiteP_z",
                 "H_HighEduP_z","H_UnempP_z","H_Income_z","H_CrowdP_z")
prognostics <- c("H_GreenPN_z","H_PublicP_z","H_ActiveP_log_z")
pca_lags    <- c("lag13_covPC1","lag13_covPC2","lag13_covPC3")
weights     <- "sw_win"
res_lags    <- c("resid_lag13")

# Build the linear‐model formula (no smooths)
f_base <- as.formula(
  paste0(
    oc, " ~ ", 
    paste(c(treat, confounders, prognostics), collapse = " + ")
  )
)
f_spac1 <- as.formula(
  paste0(
    oc, " ~ ", 
    paste(c(treat, confounders, prognostics, pca_lags), collapse = " + ")
  )
)
f_spac2 <- as.formula(
  paste0(
    oc, " ~ ", 
    paste(c(treat, confounders, prognostics, pca_lags, res_lags), collapse = " + ")
  )
)
# 3. Fit OLS (lm), *not* bam/gam ---------------------------------------------
lm_base <- lm(
  formula = f_base,
  data    = Hex_data,
  weights = Hex_data[["sw_win"]]
)
lm_spac1 <- lm(
  formula = f_spac1,
  data    = Hex_data,
  weights = Hex_data[["sw_win"]]
)
lm_spac2 <- lm(
  formula = f_spac2,
  data    = Hex_data,
  weights = Hex_data[["sw_win"]]
)

# 4. Run all LM tests --------------------------------------------------------
lm_test_base <- lm.RStests(
  model = lm_base,
  listw = list_IDW8400_13,
  test  = "all"
)
lm_test_spac1 <- lm.RStests(
  model = lm_spac1,
  listw = list_IDW8400_13,
  test  = "all"
)
lm_test_spac2 <- lm.RStests(
  model = lm_spac2,
  listw = list_IDW8400_13,
  test  = "all"
)

# 5. Inspect
print(lm_test_base)
print(lm_test_spac1)
print(lm_test_spac2)


res_base <- bind_rows(
  lapply(names(lm_test_base), function(test_name) {
    tst <- lm_test_base[[test_name]]
    data.frame(
      Test      = test_name,
      Statistic = as.numeric(tst$statistic),
      df        = as.numeric(tst$parameter),
      p_value   = as.numeric(tst$p.value),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
)

# 2. (Optional) Format p-values for readability
res_base <- res_base %>%
  mutate(
    p_value = ifelse(p_value < 2.2e-16, "<2.2e-16", format(round(p_value, 3), nsmall = 3))
  )

res_spac1 <- bind_rows(
  lapply(names(lm_test_spac1), function(test_name) {
    tst <- lm_test_spac1[[test_name]]
    data.frame(
      Test      = test_name,
      Statistic = as.numeric(tst$statistic),
      df        = as.numeric(tst$parameter),
      p_value   = as.numeric(tst$p.value),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
)

# 2. (Optional) Format p-values for readability
res_spac1 <- res_spac1 %>%
  mutate(
    p_value = ifelse(p_value < 2.2e-16, "<2.2e-16", format(round(p_value, 3), nsmall = 3))
  )


res_spac2 <- bind_rows(
  lapply(names(lm_test_spac2), function(test_name) {
    tst <- lm_test_spac2[[test_name]]
    data.frame(
      Test      = test_name,
      Statistic = as.numeric(tst$statistic),
      df        = as.numeric(tst$parameter),
      p_value   = as.numeric(tst$p.value),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
)

# 2. (Optional) Format p-values for readability
res_spac2 <- res_spac2 %>%
  mutate(
    p_value = ifelse(p_value < 2.2e-16, "<2.2e-16", format(round(p_value, 3), nsmall = 3))
  )
# 3. Inspect
print(res_base)

# 3. Save as CSV
write_csv(res_base, "Results/Lagrange/lm_test_base.csv")
write_csv(res_spac1, "Results/Lagrange/lm_test_spac1.csv")
write_csv(res_spac2, "Results/Lagrange/lm_test_spac2.csv")
