# 1. Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(car)  # needed for VIF
library(psych)
library(sf)
library(tidyverse)
library(knitr)
# Set working directory
base_dir <- "C:/Users/anapa/OneDrive - University of Cambridge/PhD - Research/P2 - Spatial Analysis/R Analysis/3rdAnalysis"
setwd(base_dir)
# Load data
Hex_lag1 <- read_rds("Data/Lag/Hex_lag1.rds")
Hex_lag13 <- read_rds("Data/Lag/Hex_lag13.rds")

# # 3. Create output folder
dir.create("Results/VIF", showWarnings = FALSE)

numeric_lag1 <- Hex_lag1 %>%
  select_if(is.numeric) %>%       # keep only numeric columns
  select(-H_ID, )
numeric_lag1 <- st_drop_geometry(numeric_lag1)

numeric_lag13 <- Hex_lag13 %>%
  select_if(is.numeric) %>%       # keep only numeric columns
  select(-H_ID, )
numeric_lag13 <- st_drop_geometry(numeric_lag13)

# 5. Custom order (adjust as needed)
id_var <- "H_ID"

treatment <- c("N1_8Area_log_z")

moderators <- c("N1_8BiodIx_z", "N1_8GreenP_rlog_z", "N1_8TreeP_z")

mediaotors <- c("H_AirPIx25_log_z")

cofounders <- c("H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z", 
                  "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z", "H_CrowdP_z")
                  
prognostics <- c("H_GreenPN_z", "H_PublicP_z", "H_ActiveP_log_z")

outcomes <- c("avgAF_z", "avgCHD_z", "avgHF_z", "avgHYP_z", "avgPAD_log_z","avgSTIA_z", 
               "avgOB18_z","avgCAN_z","avgDM17_z", "avgNDH18_z", "avgCKD18_log_z",
               "avgAST6_z","avgCOPD_z", "avgRA16_z", "avgOST50_log_z", 
               "avgDEP_z", "avgMH_log_z","avgDEM_log_z", "avgEP18_z", "avgLD_z")

cof_pca <- c("cofPC1", "cofPC2", "cofPC3", "cofPC4", "cofPC5")

cov_pca <- c("covPC1", "covPC2", "covPC3", "covPC4", "covPC5", "covPC6", "covPC7")

treatment_lag1 <- c("N1_8Area_log_z_lag1")

cofounders_lag1 <- c("H_Density_log_z_lag1", "H_Min65P_z_lag1", "H_WhiteP_z_lag1", "H_HighEduP_z_lag1",
                     "H_UnempP_z_lag1", "H_Income_z_lag1", "H_PHeatP_log_z_lag1", "H_CrowdP_z_lag1")

prognostics_lag1 <- c("H_GreenPN_z_lag1", "H_PublicP_z_lag1", "H_ActiveP_log_z_lag1")       

mediaotors_lag1 <- c("H_AirPIx25_log_z_lag1")

outcomes_lag1 <- c("avgAF_z_lag1", "avgCHD_z_lag1", "avgHF_z_lag1", "avgHYP_z_lag1", "avgPAD_log_z_lag1", "avgSTIA_z_lag1", 
                   "avgOB18_z_lag1", "avgCAN_z_lag1", "avgDM17_z_lag1", "avgNDH18_z_lag1", "avgCKD18_log_z_lag1",
                   "avgAST6_z_lag1", "avgCOPD_z_lag1", "avgRA16_z_lag1", "avgOST50_log_z_lag1", 
                   "avgDEP_z_lag1", "avgMH_log_z_lag1", "avgDEM_log_z_lag1", "avgEP18_z_lag1", "avgLD_z_lag1")

cof_pca_lag1 <- c("cofPC1_lag1", "cofPC2_lag1", "cofPC3_lag1", "cofPC4_lag1", "cofPC5_lag1")

cov_pca_lag1 <- c("covPC1_lag1", "covPC2_lag1", "covPC3_lag1", "covPC4_lag1", "covPC5_lag1", "covPC6_lag1", "covPC7_lag1")

lag1_pca_cof <- c("lag1_cofPC1", "lag1_cofPC2", "lag1_cofPC3")

lag1_pca_cov <- c("lag1_covPC1", "lag1_covPC2", "lag1_covPC3")

treatment_lag13 <- c("N1_8Area_log_z_lag13")

cofounders_lag13 <- c("H_Density_log_z_lag13", "H_Min65P_z_lag13", "H_WhiteP_z_lag13", "H_HighEduP_z_lag13",
                      "H_UnempP_z_lag13", "H_Income_z_lag13", "H_PHeatP_log_z_lag13", "H_CrowdP_z_lag13")

prognostics_lag13 <- c("H_GreenPN_z_lag13", "H_PublicP_z_lag13", "H_ActiveP_log_z_lag13")       

mediaotors_lag13 <- c("H_AirPIx25_log_z_lag13")

outcomes_lag13 <- c("avgAF_z_lag13", "avgCHD_z_lag13", "avgHF_z_lag13", "avgHYP_z_lag13", "avgPAD_log_z_lag13", "avgSTIA_z_lag13", 
                     "avgOB18_z_lag13", "avgCAN_z_lag13", "avgDM17_z_lag13", "avgNDH18_z_lag13", "avgCKD18_log_z_lag13",
                     "avgAST6_z_lag13", "avgCOPD_z_lag13", "avgRA16_z_lag13", "avgOST50_log_z_lag13", 
                     "avgDEP_z_lag13", "avgMH_log_z_lag13", "avgDEM_log_z_lag13", "avgEP18_z_lag13", "avgLD_z_lag13")

cof_pca_lag13 <- c("cofPC1_lag13", "cofPC2_lag13", "cofPC3_lag13", "cofPC4_lag13", "cofPC5_lag13")

cov_pca_lag13 <- c("covPC1_lag13", "covPC2_lag13", "covPC3_lag13", "covPC4_lag13", "covPC5_lag13", "covPC6_lag13", "covPC7_lag13")

lag13_pca_cof <- c("lag13_cofPC1", "lag13_cofPC2", "lag13_cofPC3")

lag13_pca_cov <- c("lag13_covPC1", "lag13_covPC2", "lag13_covPC3")

# 3. VIF FUNCTION GPS AREA
vif_gps_lag1 <- c(cofounders, treatment, cofounders_lag1)
vif_gps_pcalag1 <- c(cofounders, treatment, prognostics, lag1_pca_cof)
vif_gps_lag13 <- c(cofounders, treatment, cofounders_lag13)
vif_gps_pcalag13 <- c(cofounders, treatment, prognostics, lag13_pca_cov)

vif_gam_lag1 <- c(cofounders, treatment, prognostics, moderators, cofounders_lag1, prognostics_lag1)
vif_gam_pcalag1 <- c(cofounders, treatment, prognostics, moderators, lag1_pca_cov)
vif_gam_lag13 <- c(cofounders, treatment, prognostics, moderators, cofounders_lag13, prognostics_lag13)
vif_gam_pcalag13 <- c(cofounders, treatment, prognostics, moderators, lag13_pca_cov)

vif_durbin_lag1 <- c(cofounders, treatment, prognostics, treatment_lag1, cofounders_lag1, prognostics_lag1, "avgAF_z_lag1")
vif_durbin_pcalag1 <- c(cov_pca,"H_GreenPN_z", "H_GreenPN_z_lag1", treatment, treatment_lag1, cov_pca_lag1, "avgAF_z_lag1")
vif_durbin_lag13 <- c(cofounders, treatment, prognostics, treatment_lag13, cofounders_lag13, prognostics_lag13, "avgAF_z_lag13")
vif_durbin_pcalag13 <- c(cov_pca, "H_GreenPN_z", "H_GreenPN_z_lag13", treatment, treatment_lag13, cov_pca_lag13, "avgAF_z_lag13")

# VIF helper
vif_all_vars <- function(data, vars) {
  if (length(vars) < 2) stop("Need at least two variables for VIF")
  form <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse = "+")))
  lm_mod <- lm(form, data = data)
  car::vif(lm_mod)
}

vif_gps_lag1_res <- vif_all_vars(numeric_lag1, vif_gps_lag1)
print(vif_gps_lag1_res)
saveRDS(vif_gps_lag1_res, "Results/VIF/vif_gps_lag1_res.rds")

vif_gps_pcalag1_res <- vif_all_vars(numeric_lag1, vif_gps_pcalag1)
print(vif_gps_pcalag1_res)
saveRDS(vif_gps_pcalag1_res, "Results/VIF/vif_gps_pcalag1_res.rds")

vif_gam_lag1_res <- vif_all_vars(numeric_lag1, vif_gam_lag1)
print(vif_gam_lag1_res)
saveRDS(vif_gam_lag1_res, "Results/VIF/vif_gam_lag1_res.rds")

vif_gam_pcalag1_res <- vif_all_vars(numeric_lag1, vif_gam_pcalag1)
print(vif_gam_pcalag1_res)
saveRDS(vif_gam_pcalag1_res, "Results/VIF/vif_gam_pcalag1_res.rds")

vif_durbin_lag1_res <- vif_all_vars(numeric_lag1, vif_durbin_lag1)
print(vif_durbin_lag1_res)
saveRDS(vif_durbin_lag1_res, "Results/VIF/vif_durbin_lag1_res.rds")

vif_durbin_pcalag1_res <- vif_all_vars(numeric_lag1, vif_durbin_pcalag1)
print(vif_durbin_pcalag1_res)
saveRDS(vif_durbin_pcalag1_res, "Results/VIF/vif_durbin_pcalag1_res.rds")

vif_gps_lag13_res <- vif_all_vars(numeric_lag13, vif_gps_lag13)
print(vif_gps_lag13_res)
saveRDS(vif_gps_lag13_res, "Results/VIF/vif_gps_lag13_res.rds")

vif_gps_pcalag13_res <- vif_all_vars(numeric_lag13, vif_gps_pcalag13)
print(vif_gps_pcalag13_res)
saveRDS(vif_gps_pcalag13_res, "Results/VIF/vif_gps_pcalag13_res.rds")

vif_gam_lag13_res <- vif_all_vars(numeric_lag13, vif_gam_lag13)
print(vif_gam_lag13_res)
saveRDS(vif_gam_lag13_res, "Results/VIF/vif_gam_lag13_res.rds")

vif_gam_pcalag13_res <- vif_all_vars(numeric_lag13, vif_gam_pcalag13)
print(vif_gam_pcalag13_res)
saveRDS(vif_gam_pcalag13_res, "Results/VIF/vif_gam_pcalag13_res.rds")

vif_durbin_lag13_res <- vif_all_vars(numeric_lag13, vif_durbin_lag13)
print(vif_durbin_lag13_res)
saveRDS(vif_durbin_lag13_res, "Results/VIF/vif_durbin_lag13_res.rds")

vif_durbin_pcalag13_res <- vif_all_vars(numeric_lag13, vif_durbin_pcalag13)
print(vif_durbin_pcalag13_res)
saveRDS(vif_durbin_pcalag13_res, "Results/VIF/vif_durbin_pcalag13_res.rds")

###Compare results between Lags 1 and 1.3

# 1. Point to all of your VIF result files
vif_files <- tibble(
  path = list.files("Results/VIF", pattern = "^vif_.*_res\\.rds$", full.names = TRUE)
) %>%
  mutate(
    fname = basename(path),
    lag   = str_extract(fname, "lag13|lag1"),
    model = fname %>%
      str_remove("^vif_") %>%
      str_remove("lag(?:1|13)_res\\.rds$")
  )

# 2. Read each one and turn into a long tibble
vif_tbl <- vif_files %>%
  mutate(
    data = map(
      path,
      ~ read_rds(.x) %>% 
        enframe(name = "variable", value = "VIF")
    )
  ) %>%
  select(model, lag, data) %>%
  unnest(data)

# 3. (Optional) Pivot wider so you have one column per lag:
vif_wide <- vif_tbl %>%
  pivot_wider(names_from = lag, values_from = VIF) %>%
  arrange(model, desc(lag1))

# 4a. Print a kable for quick visual check
kable(vif_wide, digits = 2, caption = "Comparison of VIFs: lag1 vs lag13")

write.csv(vif_wide,  "Results/VIF/vifs_lag1_13.csv", row.names = FALSE)

