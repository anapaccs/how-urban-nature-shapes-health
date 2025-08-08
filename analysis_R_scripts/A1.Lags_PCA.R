# Load required packages
library(sf)
library(dplyr)
library(spdep)
library(tibble)
library(readr)
library(units)

# Set working directory
base_dir <- "C:/Users/anapa/OneDrive - University of Cambridge/PhD - Research/P2 - Spatial Analysis/R Analysis/3rdAnalysis"
setwd(base_dir)
# input / output
data_in  <- file.path(base_dir, "Data/Treated/HEX_data_transf.rds")
dir.create("Results/Lags", showWarnings = FALSE, recursive = TRUE)
dir.create("Data/Lag", showWarnings = FALSE, recursive = TRUE)

# 1. Load & prep ----------------------------------------------------------
HEX_sf <- read_rds(data_in)         
coords  <- st_coordinates(st_centroid(HEX_sf))


raw_vars <- c(
  # prognostics
  "H_GreenPN_z", "H_PublicP_z", "H_ActiveP_log_z",
  # cofounders
  "H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z",
  "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z", "H_CrowdP_z", 
  # treatment
  "N1_8Area_log_z",
  #moderators
  "N1_8BiodIx_z", "N1_8GreenP_rlog_z", "N1_8TreeP_z", 
  #outcomes
  "avgAF_z", "avgCHD_z", "avgHF_z", "avgHYP_z", "avgPAD_log_z","avgSTIA_z", 
  "avgOB18_z","avgCAN_z","avgDM17_z", "avgNDH18_z", "avgCKD18_log_z",
  "avgAST6_z","avgCOPD_z", "avgRA16_z", "avgOST50_log_z", 
  "avgDEP_z", "avgMH_log_z","avgDEM_log_z", "avgEP18_z", "avgLD_z"
  )

HEX_dt <- HEX_sf %>%
  filter(!is.na(H_Density_log_z)) %>%
  select(H_ID, all_of(raw_vars))

n <- nrow(HEX_dt)

# 2. PCA on raw variables ------------------------------------------------------
#2.1 Cofounders
cof_vars <- c(
  "H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z", 
  "H_UnempP_z", "H_Income_z",  "H_PHeatP_log_z", "H_CrowdP_z")

cc <- complete.cases(HEX_dt[, cof_vars])
pca_raw <- prcomp(HEX_dt[cc, cof_vars], center = TRUE, scale. = TRUE)
cumvar  <- cumsum(pca_raw$sdev^2 / sum(pca_raw$sdev^2))
k_raw   <- which(cumvar >= 0.90)[1]

# scores (with NAs for incomplete rows)
rawPC <- matrix(NA, nrow = n, ncol = k_raw)
rawPC[cc, ] <- pca_raw$x[, 1:k_raw]
colnames(rawPC) <- paste0("cofPC", 1:k_raw)

HEX_dt <- bind_cols(HEX_dt, as_tibble(rawPC))

#2.1 Covariates
cov_vars <- c(
  # confounders + prognostics
  "H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z", 
  "H_UnempP_z",      "H_Income_z",  "H_PHeatP_log_z", "H_CrowdP_z",
  "H_PublicP_z",     "H_ActiveP_log_z")

cc_cov <- complete.cases(HEX_dt[, cov_vars])
pca_cov <- prcomp(HEX_dt[cc_cov, cov_vars], center = TRUE, scale. = TRUE)
cumcov  <- cumsum(pca_cov$sdev^2 / sum(pca_cov$sdev^2))
k_cov   <- which(cumcov >= 0.90)[1]

# scores (with NAs for incomplete rows)
covPC <- matrix(NA, nrow = n, ncol = k_cov)
covPC[cc, ] <- pca_cov$x[, 1:k_cov]
colnames(covPC) <- paste0("covPC", 1:k_cov)

HEX_dt <- bind_cols(HEX_dt, as_tibble(covPC))

raw_PCs <- c("cofPC1", "cofPC2", "cofPC3", "cofPC4", "cofPC5", "covPC1",
                 "covPC2", "covPC3", "covPC4", "covPC5", "covPC6", "covPC7")  

# 3 build one dataset per weighting --------------------------------------------

# 3.1 On the FULL sf (n≈15 000), compute 1/d & 1/d^1.3 lags
# compute centroids & pairwise distances
cent_full <- st_centroid(HEX_sf)
d_full    <- st_distance(cent_full) %>% drop_units()

# zero out self‐distances 
diag(d_full) <- NA

# 1/d weights, row‐standardized
w1_full  <- d_full^-1;    w1_full[!is.finite(w1_full)] <- 0;    w1_full <- w1_full / rowSums(w1_full)
# 1/d^1.3 weights, row‐standardized
w13_full <- d_full^-1.3;  w13_full[!is.finite(w13_full)] <- 0;  w13_full <- w13_full / rowSums(w13_full)

# attach lag columns back onto HEX_sf
HEX_sf <- HEX_sf %>% mutate(
  H_GreenPN_z_lag1    = as.numeric(w1_full  %*% H_GreenPN_z),
  H_GreenPN_z_lag13  = as.numeric(w13_full %*% H_GreenPN_z),
  H_AirPIx25_log_z_lag1   = as.numeric(w1_full  %*% H_AirPIx25_log_z),
  H_AirPIx25_log_z_lag13 = as.numeric(w13_full %*% H_AirPIx25_log_z)
)

# 3.2 Filter down to 12 469 rows 
HEX_filt <- HEX_sf %>% 
  filter(!is.na(H_Density_log_z))

# Build df for the rest of the workflow
df_base <- HEX_filt %>%
  st_drop_geometry() %>% 
  select(H_ID, all_of(raw_vars)) %>%                 
  left_join(
    HEX_dt %>% select(H_ID, all_of(raw_PCs)), 
    by = "H_ID"
  ) %>%
  left_join(
    HEX_sf %>% st_drop_geometry() %>% 
      select(H_ID, 
             H_GreenPN_z_lag1, H_GreenPN_z_lag13,
             H_AirPIx25_log_z_lag1, H_AirPIx25_log_z_lag13),
    by = "H_ID"
  )

# 3.4 Lag all other vars on df_base
# Get the filtered centroids & full distance matrix
cent_filt <- st_centroid(HEX_filt)
d_inf     <- st_distance(cent_filt) %>% drop_units()
diag(d_inf) <- NA          # zero out self‐distances

# Make 1/d and 1/d^1.3, row‐standardized
w1_filt   <- d_inf^-1
w1_filt[!is.finite(w1_filt)] <- 0
w1_filt   <- w1_filt / rowSums(w1_filt)

w13_filt  <- d_inf^-1.3
w13_filt[!is.finite(w13_filt)] <- 0
w13_filt  <- w13_filt / rowSums(w13_filt)

# Pick vars to lag
vars_to_lag <- c("H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z", "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z",
                 "H_CrowdP_z", "H_PublicP_z", "H_ActiveP_log_z", "N1_8Area_log_z", "N1_8BiodIx_z", "N1_8GreenP_rlog_z", 
                 "N1_8TreeP_z", "avgAF_z", "avgCHD_z", "avgHF_z", "avgHYP_z", "avgPAD_log_z", "avgSTIA_z", "avgOB18_z", "avgCAN_z", 
                 "avgDM17_z", "avgNDH18_z", "avgCKD18_log_z", "avgAST6_z", "avgCOPD_z", "avgRA16_z", "avgOST50_log_z", "avgDEP_z", 
                 "avgMH_log_z", "avgDEM_log_z", "avgEP18_z", "avgLD_z", "cofPC1", "cofPC2", "cofPC3", "cofPC4", "cofPC5", "covPC1",
                 "covPC2", "covPC3", "covPC4", "covPC5", "covPC6", "covPC7")  

# Matrix‐multiply to get lag1 and lag1.3 
X <- as.matrix(df_base[, vars_to_lag])

# Compute all the 1/d lags at once (n × n) %*% (n × p) → (n × p)
lag1_mat  <- w1_filt  %*% X
lag13_mat <- w13_filt %*% X

colnames(lag1_mat)  <- paste0(vars_to_lag, "_lag1")
colnames(lag13_mat) <- paste0(vars_to_lag, "_lag13")

df_lagged <- bind_cols(
  df_base,
  as_tibble(lag1_mat),
  as_tibble(lag13_mat)
)

colnames(df_lagged)

#----------------------------------
# Step 4: PCA on Lagged Raw Variables (e.g. from IDW13)
#----------------------------------
cof_vars_lag1 <- c(
  "H_Density_log_z_lag1", "H_Min65P_z_lag1", "H_WhiteP_z_lag1", "H_HighEduP_z_lag1", 
  "H_UnempP_z_lag1", "H_Income_z_lag1",  "H_PHeatP_log_z_lag1", "H_CrowdP_z_lag1")

cc1         <- complete.cases(df_lagged[, cof_vars_lag1])
pca_lag1    <- prcomp(df_lagged[cc1, cof_vars_lag1], center = TRUE, scale. = TRUE)
varcum1     <- cumsum(pca_lag1$sdev^2 / sum(pca_lag1$sdev^2))
k_lag1      <- which(varcum1 >= 0.90)[1]

lagPC1      <- matrix(NA, nrow = nrow(df_lagged), ncol = k_lag1)
lagPC1[cc1, ] <- pca_lag1$x[, 1:k_lag1]
colnames(lagPC1) <- paste0("lag1_cofPC", 1:k_lag1)

df_lagged <- bind_cols(df_lagged, as_tibble(lagPC1))

cof_vars_lag13 <- c(
  "H_Density_log_z_lag13", "H_Min65P_z_lag13", "H_WhiteP_z_lag13", "H_HighEduP_z_lag13", 
  "H_UnempP_z_lag13", "H_Income_z_lag13",  "H_PHeatP_log_z_lag13", "H_CrowdP_z_lag13")

cc2         <- complete.cases(df_lagged[, cof_vars_lag13])
pca_lag2    <- prcomp(df_lagged[cc2, cof_vars_lag13], center = TRUE, scale. = TRUE)
varcum2     <- cumsum(pca_lag2$sdev^2 / sum(pca_lag2$sdev^2))
k_lag2      <- which(varcum2 >= 0.90)[1]

lagPC2      <- matrix(NA, nrow = nrow(df_lagged), ncol = k_lag2)
lagPC2[cc2, ] <- pca_lag2$x[, 1:k_lag2]
colnames(lagPC2) <- paste0("lag13_cofPC", 1:k_lag2)

df_lagged <- bind_cols(df_lagged, as_tibble(lagPC2))

cov_vars_lag1 <- c(
  "H_Density_log_z_lag1", "H_Min65P_z_lag1", "H_WhiteP_z_lag1", "H_HighEduP_z_lag1", 
  "H_UnempP_z_lag1", "H_Income_z_lag1",  "H_PHeatP_log_z_lag1", "H_CrowdP_z_lag1")

cc3         <- complete.cases(df_lagged[, cov_vars_lag1])
pca_lag3    <- prcomp(df_lagged[cc3, cov_vars_lag1], center = TRUE, scale. = TRUE)
varcum3     <- cumsum(pca_lag3$sdev^2 / sum(pca_lag3$sdev^2))
k_lag3      <- which(varcum3 >= 0.90)[1]

lagPC3      <- matrix(NA, nrow = nrow(df_lagged), ncol = k_lag3)
lagPC3[cc3, ] <- pca_lag3$x[, 1:k_lag3]
colnames(lagPC3) <- paste0("lag1_covPC", 1:k_lag3)

df_lagged <- bind_cols(df_lagged, as_tibble(lagPC3))

cov_vars_lag13 <- c(
  "H_Density_log_z_lag13", "H_Min65P_z_lag13", "H_WhiteP_z_lag13", "H_HighEduP_z_lag13", 
  "H_UnempP_z_lag13", "H_Income_z_lag13",  "H_PHeatP_log_z_lag13", "H_CrowdP_z_lag13")

cc4         <- complete.cases(df_lagged[, cov_vars_lag13])
pca_lag4    <- prcomp(df_lagged[cc4, cov_vars_lag13], center = TRUE, scale. = TRUE)
varcum4     <- cumsum(pca_lag4$sdev^2 / sum(pca_lag4$sdev^2))
k_lag4      <- which(varcum4 >= 0.90)[1]

lagPC4      <- matrix(NA, nrow = nrow(df_lagged), ncol = k_lag4)
lagPC4[cc1, ] <- pca_lag4$x[, 1:k_lag4]
colnames(lagPC4) <- paste0("lag13_covPC", 1:k_lag4)

df_lagged <- bind_cols(df_lagged, as_tibble(lagPC4))

#----------------------------------
# Step 5: Save Results
#----------------------------------
colnames(df_lagged)

Hex_lag1 <- df_lagged %>%
  select("H_ID",
         # treatment
         "N1_8Area_log_z",
         #moderators
         "N1_8BiodIx_z", "N1_8GreenP_rlog_z", "N1_8TreeP_z", 
         # cofounders
         "H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z",
         "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z", "H_CrowdP_z", 
         # prognostics
         "H_GreenPN_z", "H_PublicP_z", "H_ActiveP_log_z",
         # mediator
         "H_AirPIx25_log_z",
         #outcomes
         "avgAF_z", "avgCHD_z", "avgHF_z", "avgHYP_z", "avgPAD_log_z","avgSTIA_z", 
         "avgOB18_z","avgCAN_z","avgDM17_z", "avgNDH18_z", "avgCKD18_log_z",
         "avgAST6_z","avgCOPD_z", "avgRA16_z", "avgOST50_log_z", 
         "avgDEP_z", "avgMH_log_z","avgDEM_log_z", "avgEP18_z", "avgLD_z",
         #raw_pca_cof
         "cofPC1", "cofPC2", "cofPC3", "cofPC4", "cofPC5", 
         #raw_pca_cov
         "covPC1", "covPC2", "covPC3", "covPC4", "covPC5", "covPC6", "covPC7",
         # treatmentlag1
         "N1_8Area_log_z_lag1",
         #moderatorslag1
         "N1_8BiodIx_z_lag1", "N1_8GreenP_rlog_z_lag1", "N1_8TreeP_z_lag1", 
         # cofounderslag1
         "H_Density_log_z_lag1", "H_Min65P_z_lag1", "H_WhiteP_z_lag1", "H_HighEduP_z_lag1",
         "H_UnempP_z_lag1", "H_Income_z_lag1", "H_PHeatP_log_z_lag1", "H_CrowdP_z_lag1", 
         # prognosticslag1
         "H_GreenPN_z_lag1", "H_PublicP_z_lag1", "H_ActiveP_log_z_lag1",
         # mediatorlag1
         "H_AirPIx25_log_z_lag1",
         #outcomeslag1
         "avgAF_z_lag1", "avgCHD_z_lag1", "avgHF_z_lag1", "avgHYP_z_lag1", "avgPAD_log_z_lag1", "avgSTIA_z_lag1", 
         "avgOB18_z_lag1", "avgCAN_z_lag1", "avgDM17_z_lag1", "avgNDH18_z_lag1", "avgCKD18_log_z_lag1",
         "avgAST6_z_lag1", "avgCOPD_z_lag1", "avgRA16_z_lag1", "avgOST50_log_z_lag1", 
         "avgDEP_z_lag1", "avgMH_log_z_lag1", "avgDEM_log_z_lag1", "avgEP18_z_lag1", "avgLD_z_lag1",
         #raw_pca_coflag1
         "cofPC1_lag1", "cofPC2_lag1", "cofPC3_lag1", "cofPC4_lag1", "cofPC5_lag1", 
         #raw_pca_covlag1
         "covPC1_lag1", "covPC2_lag1", "covPC3_lag1", "covPC4_lag1", "covPC5_lag1", "covPC6_lag1", "covPC7_lag1",
         #pca_lag1_cof
         "lag1_cofPC1", "lag1_cofPC2", "lag1_cofPC3",
         #pca_lag1_covlag
         "lag1_covPC1", "lag1_covPC2", "lag1_covPC3")

Hex_lag13 <- df_lagged %>%
  select("H_ID",
         # treatment
         "N1_8Area_log_z",
         #moderators
         "N1_8BiodIx_z", "N1_8GreenP_rlog_z", "N1_8TreeP_z", 
         # cofounders
         "H_Density_log_z", "H_Min65P_z", "H_WhiteP_z", "H_HighEduP_z",
         "H_UnempP_z", "H_Income_z", "H_PHeatP_log_z", "H_CrowdP_z", 
         # prognostics
         "H_GreenPN_z", "H_PublicP_z", "H_ActiveP_log_z",
         # mediator
         "H_AirPIx25_log_z",
         #outcomes
         "avgAF_z", "avgCHD_z", "avgHF_z", "avgHYP_z", "avgPAD_log_z","avgSTIA_z", 
         "avgOB18_z","avgCAN_z","avgDM17_z", "avgNDH18_z", "avgCKD18_log_z",
         "avgAST6_z","avgCOPD_z", "avgRA16_z", "avgOST50_log_z", 
         "avgDEP_z", "avgMH_log_z","avgDEM_log_z", "avgEP18_z", "avgLD_z",
         #raw_pca_cof
         "cofPC1", "cofPC2", "cofPC3", "cofPC4", "cofPC5", 
         #raw_pca_cov
         "covPC1", "covPC2", "covPC3", "covPC4", "covPC5", "covPC6", "covPC7",
         # treatmentlag13
         "N1_8Area_log_z_lag13",
         #moderatorslag13
         "N1_8BiodIx_z_lag13", "N1_8GreenP_rlog_z_lag13", "N1_8TreeP_z_lag13", 
         # cofounderslag13
         "H_Density_log_z_lag13", "H_Min65P_z_lag13", "H_WhiteP_z_lag13", "H_HighEduP_z_lag13",
         "H_UnempP_z_lag13", "H_Income_z_lag13", "H_PHeatP_log_z_lag13", "H_CrowdP_z_lag13", 
         # prognosticslag13
         "H_GreenPN_z_lag13", "H_PublicP_z_lag13", "H_ActiveP_log_z_lag13",
         # mediatorlag13
         "H_AirPIx25_log_z_lag13",
         #outcomeslag13
         "avgAF_z_lag13", "avgCHD_z_lag13", "avgHF_z_lag13", "avgHYP_z_lag13", "avgPAD_log_z_lag13", "avgSTIA_z_lag13", 
         "avgOB18_z_lag13", "avgCAN_z_lag13", "avgDM17_z_lag13", "avgNDH18_z_lag13", "avgCKD18_log_z_lag13",
         "avgAST6_z_lag13", "avgCOPD_z_lag13", "avgRA16_z_lag13", "avgOST50_log_z_lag13", 
         "avgDEP_z_lag13", "avgMH_log_z_lag13", "avgDEM_log_z_lag13", "avgEP18_z_lag13", "avgLD_z_lag13",
         #raw_pca_coflag13
         "cofPC1_lag13", "cofPC2_lag13", "cofPC3_lag13", "cofPC4_lag13", "cofPC5_lag13", 
         #raw_pca_covlag13
         "covPC1_lag13", "covPC2_lag13", "covPC3_lag13", "covPC4_lag13", "covPC5_lag13", "covPC6_lag13", "covPC7_lag13",
         #pca_lag13_cof
         "lag13_cofPC1", "lag13_cofPC2", "lag13_cofPC3",
         #pca_lag13_covlag
         "lag13_covPC1", "lag13_covPC2", "lag13_covPC3")



# --- g) Save your two “infinite‐range” datasets
write_rds(Hex_lag1,  "Data/Lag/Hex_lag1.rds")
write.csv(Hex_lag1,  "Data/Lag/Hex_lag1.csv", row.names = FALSE)
write_rds(Hex_lag13, "Data/Lag/Hex_lag13.rds")
write.csv(Hex_lag13,  "Data/Lag/Hex_lag13.csv", row.names = FALSE)

#----------------------------------
# Step 2: Create Spatial Weight Matrices
#----------------------------------
# Assuming you already have coords and idw:
coords <- st_coordinates(st_centroid(HEX_dt))  

# 2a. Distance-band matrices
idw2800 <- dnearneigh(coords, 0, 2800)
idw5600 <- dnearneigh(coords, 0, 5600)
idw8400 <- dnearneigh(coords, 0, 8400)

#distances <- nbdists(idw, coords)
distances2800 <- nbdists(idw2800, coords)
distances5600 <- nbdists(idw5600, coords)
distances8400 <- nbdists(idw8400, coords)

idw2800_w <- nb2listw(idw2800, glist = lapply(distances2800, function(d) 1 / d), style = "W", zero.policy = TRUE)
idw5600_w <- nb2listw(idw5600, glist = lapply(distances5600, function(d) 1 / d), style = "W", zero.policy = TRUE)
idw8400_w <- nb2listw(idw8400, glist = lapply(distances8400, function(d) 1 / d), style = "W", zero.policy = TRUE)

idw2800_13w <- nb2listw(idw2800, glist = lapply(distances2800, function(d) 1 / d^1.3), style = "W", zero.policy = TRUE)
idw5600_13w <- nb2listw(idw5600, glist = lapply(distances5600, function(d) 1 / d^1.3), style = "W", zero.policy = TRUE)
idw8400_13w <- nb2listw(idw8400, glist = lapply(distances8400, function(d) 1 / d^1.3), style = "W", zero.policy = TRUE)

# Save weight matrices
saveRDS(idw2800_w, "Results/Lags/idw2800_w.rds")
saveRDS(idw2800_13w, "Results/Lags/idw2800_13w.rds")
saveRDS(idw5600_w, "Results/Lags/idw5600_w.rds")
saveRDS(idw5600_13w, "Results/Lags/idw5600_13w.rds")
saveRDS(idw8400_w, "Results/Lags/idw8400_w.rds")
saveRDS(idw8400_13w, "Results/Lags/idw8400_13w.rds")

saveRDS(idw2800_13w, "Results/Lags/idw2800_13w_filt.rds")
saveRDS(idw5600_13w, "Results/Lags/idw5600_13w_filt.rds")
saveRDS(idw8400_13w, "Results/Lags/idw8400_13w_filt.rds")
