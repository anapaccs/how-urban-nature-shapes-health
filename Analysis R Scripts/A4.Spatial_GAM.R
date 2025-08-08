# ------------------------------------------------------------------------
# Full GAM pipeline: setup, fit models, plots, maps, summaries
# ------------------------------------------------------------------------

# 0) LIBRARIES & WORKING DIRECTORY ---------------------------------------
library(tidyverse)   
library(mgcv)         
library(broom)       
library(sf)          
library(gratia) 
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)  
library(units)
library(gridExtra)
library(RColorBrewer)  
library(paletteer)

# Use in a ggplot2 chart:
scale_colour_paletteer_d("LaCroixColoR::Berry")
scale_fill_paletteer_d("LaCroixColoR::Berry")

base_dir <- "C:/Users/anapa/OneDrive - University of Cambridge/PhD - Research/P2 - Spatial Analysis/R Analysis/3rdAnalysis"
setwd(base_dir)

# create output folders
dir.create("Results/Gamfs",         recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Gamfs/Plots",   recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Gamfs/Maps",    recursive = TRUE, showWarnings = FALSE)

# 1) LOAD & ENSURE sf GEOMETRY -------------------------------------------
Hex_data     <- read_rds("Data/Mod/Hex_data.rds")
hex_sp       <- st_read("Data/Treated/HEX_data_transf.gpkg", quiet=TRUE)
centroids <- st_centroid(hex_sp)
centroids_df <- centroids %>%
  mutate(
    lon = st_coordinates(centroids)[,1],
    lat = st_coordinates(centroids)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(H_ID, lon, lat)

HEX_new  <- Hex_data %>%
  left_join(centroids_df, by = "H_ID")

hex_coords <- hex_sp %>%
  st_centroid() %>%
  transmute(H_ID,
            lon = st_coordinates(.)[,1],
            lat = st_coordinates(.)[,2])

Hex_data <- Hex_data %>%
  left_join(hex_coords, by="H_ID") 

HEX_sf <- Hex_data

# drop geometry for modeling
HEX_df <- st_drop_geometry(HEX_sf)
HEX_df <- HEX_df %>%
  # 1) create a single lon/lat using coalesce()
  mutate(
    lon = coalesce(lon.x, lon.y),
    lat = coalesce(lat.x, lat.y)
  ) %>%
  # 2) drop all the old, now-redundant columns
  select(
    -lon.x, -lat.x,
    -lon.y, -lat.y,
    -geom.x, -geom.y  # if these were also spurious
  )

# 2) PRECOMPUTE SPATIAL-WEIGHTS W13 ---------------------------------------
HEX_sf <- st_as_sf(Hex_data)
cent   <- st_centroid(HEX_sf)
d_inf  <- st_distance(cent) %>% drop_units()
diag(d_inf) <- NA
W13    <- d_inf^(-1.3); W13[!is.finite(W13)] <- 0; W13 <- W13/rowSums(W13)

# 3) SETUP: outcomes, transforms, covariates -----------------------------
outcomes_df <- tibble(
  outcome   = c(
    "avgAF_z","avgCHD_z","avgHF_z","avgHYP_z","avgPAD_log_z",
    "avgSTIA_z","avgOB18_z","avgCAN_z","avgDM17_z","avgNDH18_z",
    "avgCKD18_log_z","avgAST6_z","avgCOPD_z","avgRA16_z",
    "avgOST50_log_z","avgDEP_z","avgMH_log_z","avgDEM_log_z",
    "avgEP18_z","avgLD_z"
  ),
  raw_name  = c(
    "avgAF","avgCHD","avgHF","avgHYP","avgPAD","avgSTIA","avgOB18","avgCAN",
    "avgDM17","avgNDH18","avgCKD18","avgAST6","avgCOPD","avgRA16","avgOST50",
    "avgDEP","avgMH","avgDEM","avgEP18","avgLD"
  ),
  transform = if_else(str_detect(outcome, "_log_z"), "log1p", "identity")
)

# compute μ/σ on transformed scale for each outcome
outcomes_stats <- outcomes_df %>%
  rowwise() %>%
  mutate(
    μ_y = if (transform=="log1p") mean(log1p(HEX_df[[raw_name]]), na.rm=TRUE)
    else                   mean(      HEX_df[[raw_name]], na.rm=TRUE),
    σ_y = if (transform=="log1p") sd(  log1p(HEX_df[[raw_name]]), na.rm=TRUE)
    else                   sd(        HEX_df[[raw_name]], na.rm=TRUE)
  ) %>% ungroup()

# treatment variable & its μ/σ on log1p scale
treat       <- "N1_8Area_log_z"
treat_stats <- tibble(
  μ_t = mean(log1p(HEX_df$N1_8Area), na.rm=TRUE),
  σ_t = sd(  log1p(HEX_df$N1_8Area), na.rm=TRUE)
)

# inverse‐transform functions
inv_id       <- function(z, μ, σ)       z * σ + μ
inv_log1p    <- function(z, μ, σ) exp(z * σ + μ) - 1

# covariates
confounders <- c("H_Density_log_z","H_Min65P_z","H_WhiteP_z",
                 "H_HighEduP_z","H_UnempP_z","H_Income_z","H_CrowdP_z")
prognostics <- c("H_GreenPN_z","H_PublicP_z","H_ActiveP_log_z")
pca_lags    <- c("lag13_covPC1","lag13_covPC2","lag13_covPC3")
moderators  <- c("N1_8BiodIx_z","N1_8GreenP_rlog_z","N1_8TreeP_z")
weights     <- "sw_win"

cov_means <- Hex_data %>%
  summarise(across(all_of(c(confounders, prognostics, pca_lags)),
                   ~ mean(.x, na.rm = TRUE)))

# grid of 200 z-values for dose–response
t_grid <- seq(min(HEX_df[[treat]], na.rm=TRUE),
              max(HEX_df[[treat]], na.rm=TRUE),
              length.out = 200)

# “other covariates” held at their means
om      <- HEX_df %>%
  summarise(across(all_of(c(confounders, prognostics, pca_lags)),
                   mean, na.rm=TRUE))

# containers
all_summaries        <- tibble()
dose_z_plots         <- list()
dose_orig_plots      <- list()
dose_zoom_plots      <- list()
int1_plots           <- list()
int2_plots           <- list()
int3_plots           <- list()

# helper to summarise GAM
summarise_gam <- function(mod, nm) {
  sm <- summary(mod)
  
  # parametric terms
  ptab <- as_tibble(sm$p.table, rownames = "term") %>%
    rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value   = `Pr(>|t|)`
    ) %>%
    mutate(type = "parametric")
  
  # smooth terms
  stab <- as_tibble(sm$s.table, rownames = "term") %>%
    rename(
      edf       = edf,
      ref.df    = `Ref.df`,
      statistic = `F`,
      p.value   = `p-value`
    ) %>%
    mutate(type = "smooth")
  
  # basic model metrics from the summary
  metrics_basic <- tibble(
    outcome  = nm,
    r.sq     = sm$r.sq,
    dev.expl = sm$dev.expl,
    n.obs    = length(mod$y)
  )
  
  # info criteria & smoothing criterion from summary
  diagnostics <- tibble(
    AIC      = AIC(mod),
    GCV_UBRE = sm$gcv.ubre,
    edf_total = sum(sm$edf) + nrow(sm$p.table),
    max_concurv = {
      conc_mat <- concurvity(mod, full = FALSE)$estimate
      if (length(conc_mat) > 1) max(conc_mat[lower.tri(conc_mat)], na.rm = TRUE) else NA
    }
  )
  
  bind_rows(ptab, stab) %>%
    mutate(outcome = nm) %>%
    left_join(metrics_basic, by = "outcome") %>%
    bind_cols(diagnostics)
}

# 4) LOOP OVER OUTCOMES ---------------------------------------------------
for (i in seq_len(nrow(outcomes_stats))) {
  oc   <- outcomes_stats$outcome[i]
  raw  <- outcomes_stats$raw_name[i]
  trfm <- outcomes_stats$transform[i]
  μ_y    <- outcomes_stats$μ_y[i]
  σ_y <- outcomes_stats$σ_y[i]
  
  # 1) Filter to complete cases
  req    <- c(oc, treat, confounders, prognostics, pca_lags)
  keep   <- complete.cases(HEX_df[, req])
  HEX_df_sub <- HEX_df[keep, ]
  HEX_sf_sub <- HEX_sf[keep, ]
  
  message("▶️  ", oc, " (n=", nrow(HEX_df_sub), " rows)")
  
  # 2) Fit the base GAM on the filtered data
  f_base   <- as.formula(
    paste0(oc, " ~ s(", treat, ", k=12, bs='ts') + ",
           paste(c(confounders, prognostics, pca_lags), collapse = " + "))
  )
  mod_base <- bam(
    f_base,
    data    = HEX_df_sub,
    weights = HEX_df_sub[[weights]],
    method  = "REML", gamma = 1.4,
    select  = TRUE, discrete = TRUE
  )
  saveRDS(mod_base, file = paste0("Results/Gamfs/mod_base_", oc, ".rds"))
  
  # 4) Compute the lag of the base residuals and add to the sub‐DF
  resid_z   <- residuals(mod_base, type = "response")
  W13_sub   <- W13[keep, keep]
  HEX_df_sub$resid_lag13 <- as.numeric(W13_sub %*% resid_z)
  summary(HEX_df_sub$resid_lag13)
  
  # 5) Fit the spatial GAM on the *same* filtered data
  f_spat   <- update(f_base, ". ~ . + resid_lag13")
  mod_spat <- bam(
    f_spat,
    data    = HEX_df_sub,
    weights = HEX_df_sub[[weights]],
    method  = "REML", gamma = 1.4,
    discrete= TRUE
  )
  
  saveRDS(mod_spat, file = paste0("Results/Gamfs/mod_spat_", oc, ".rds"))
  
  # Summaries
  sb <- summarise_gam(mod_base, oc)  %>% mutate(model="base")
  ss <- summarise_gam(mod_spat, oc)  %>% mutate(model="spatial")
  all_summaries <- bind_rows(all_summaries, sb, ss)
  
  if(oc == "avgAF_z") {
    print(summary(mod_spat)$p.table["resid_lag13", , drop=FALSE])
  }
  
  # Diagnostics
  png(paste0("Results/Gamfs/Plots/", oc, "_gamcheck_spat.png"),
      width = 800, height = 600, res = 100)
  gam.check(mod_spat, rep = 0)
  dev.off()
  
  png(paste0("Results/Gamfs/Plots/", oc, "_gamcheck_base.png"),
      width = 800, height = 600, res = 100)
  gam.check(mod_base, rep = 0)
  dev.off()
  
  # Prepare grid & predict
  base_grid_z <- tibble(t_std = t_grid) %>%
    bind_cols(om[rep(1, length(t_grid)), ]) %>%
    mutate(resid_lag13 = 0)
  newdata_z <- rename(base_grid_z, !!treat := t_std)
  
  pr <- predict(mod_spat, newdata=newdata_z, type="response", se.fit=TRUE)
  base_grid_z <- base_grid_z %>%
    mutate(pred_z = pr$fit, se_z = pr$se.fit,
           ci_low = pred_z - 2*se_z,
           ci_high= pred_z + 2*se_z)
  
  # Dose–response (z)
  p_z <- ggplot(base_grid_z, aes(x=t_std)) +
    geom_ribbon(aes(ymin=ci_low, ymax=ci_high), fill="grey70", alpha=0.4) +
    geom_line(aes(y=pred_z), linewidth = 1) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    labs(title=oc, subtitle="dose–response (z)",
         x=treat, y=paste0(oc," (z)")) +
    theme_minimal()
  
  dose_z_plots[[oc]] <- p_z
  ggsave(paste0("Results/Gamfs/Plots/", oc, "_dose_z_sp.png"),
         plot=p_z, width=6, height=4, dpi=300)
  
  # Dose–response (original)
  base_grid_orig <- base_grid_z %>%
    mutate(
      t_orig   = inv_log1p(t_std, treat_stats$μ_t, treat_stats$σ_t),
      pred_raw = if (trfm=="log1p")
        inv_log1p(pred_z, μ_y, σ_y)
      else
        inv_id(pred_z, μ_y, σ_y)
    )
  p_orig <- ggplot(base_grid_orig, aes(x=t_orig, y=pred_raw)) +
    geom_line(linewidth = 1) +
    labs(x="Distance‐decayed Area (ha)", y=raw,
         title=paste0(oc, ": dose–response (orig)")) +
    theme_minimal()
  dose_orig_plots[[oc]] <- p_orig
  ggsave(paste0("Results/Gamfs/Plots/", oc, "_dose_orig_sp.png"),
         plot=p_orig, width=6, height=4, dpi=300)
  
  # Zoomed original
  g_zoom <- p_orig +
    coord_cartesian(xlim=c(0,250)) +
    labs(subtitle="x ∈ [0, 250]")
  dose_zoom_plots[[oc]] <- g_zoom
  ggsave(paste0("Results/Gamfs/Plots/", oc, "_dose_orig_sp_0-250.png"),
         plot=g_zoom, width=6, height=4, dpi=300)
  
  # Moderation interactions (z‐scale)
  for (j in seq_along(moderators)) {
    mod_var <- moderators[j]
    # build the interaction formula properly
    f_int <- update(
      f_spat,
      as.formula(paste0(". ~ . + ", mod_var, " + ", treat, ":", mod_var))
    )
    mod_i <- bam(f_int,
                 data     = HEX_df_sub,
                 weights  = HEX_df_sub[[weights]],
                 method   = "REML",
                 discrete = TRUE)
    saveRDS(mod_i,
            file = paste0("Results/Gamfs/mod_spat_int_",
                          oc, "_", mod_var, ".rds"))
    
    qs <- quantile(HEX_df[[mod_var]], c(.1, .5, .9), na.rm=TRUE)
    grid_i <- expand_grid(t_std = t_grid, m_std = qs) %>%
      bind_cols(om[rep(1, length(t_grid)*3), ]) %>%
      mutate(resid_lag13 = 0) %>%
      mutate(pred_z = predict(mod_i,
                              newdata = rename(.,
                                               !!treat   := t_std,
                                               !!mod_var := m_std),
                              type    = "response"))
    p_i <- ggplot(grid_i,
                  aes(x = t_std, y = pred_z, color = factor(m_std))) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      scale_color_discrete(name   = mod_var,
                           labels = c("10th pct","50th pct","90th pct")) +
      labs(x     = paste0(treat, " (z)"),
           y     = paste0(oc,    " (z)"),
           title = paste0(oc,    " × ", mod_var)) +
      theme_minimal()
    if (j==1) int1_plots[[oc]] <- p_i
    if (j==2) int2_plots[[oc]] <- p_i
    if (j==3) int3_plots[[oc]] <- p_i
    ggsave(paste0("Results/Gamfs/Plots/", oc, "_int_",
                  mod_var, "_sp.png"),
           plot   = p_i, width = 6, height = 4, dpi = 300)
  }
  
  
  message("✔ Done: ", oc, "\n")
}

# 5) AGGREGATE & EXPORT SUMMARIES -----------------------------------------
sum_df <- bind_rows(all_summaries)
write_csv(sum_df, "Results/Gamfs/models_summaries_long.csv")

# pivot to wide
# keep only the spatial‐model rows
sum_df <- sum_df %>% 
  filter(model == "spatial")
sum_spat_unique <- sum_df %>%
  distinct(outcome, term, .keep_all = TRUE)

param_w <- sum_spat_unique %>% 
  filter(type=="parametric") %>%
  pivot_wider(
    id_cols     = outcome,
    names_from  = term,
    values_from = c(estimate, std.error, statistic, p.value),
    names_glue  = "{term}_{.value}"
  )

smooth_w <- sum_spat_unique %>% 
  filter(type=="smooth") %>%
  pivot_wider(
    id_cols     = outcome,
    names_from  = term,
    values_from = c(edf, `ref.df`, statistic, p.value),
    names_glue  = "{term}_{.value}"
  )

# metrics
metrics_w <- sum_spat_unique %>% 
  distinct(outcome, r.sq, dev.expl, n.obs, AIC, edf_total, max_concurv)

# join together
wide_df <- 
  metrics_w %>% 
  left_join(param_w,  by = "outcome") %>% 
  left_join(smooth_w, by = "outcome")

write_csv(wide_df, "Results/Gamfs/model_wide_spatial.csv")

# keep only the spatial‐model rows
sum_base <- sum_df %>% 
  filter(model == "base")
sum_base_unique <- sum_base  %>%
  distinct(outcome, term, .keep_all = TRUE)

param_b <- sum_base_unique %>% 
  filter(type=="parametric") %>%
  pivot_wider(
    id_cols     = outcome,
    names_from  = term,
    values_from = c(estimate, std.error, statistic, p.value),
    names_glue  = "{term}_{.value}"
  )

smooth_b <- sum_base_unique %>% 
  filter(type=="smooth") %>%
  pivot_wider(
    id_cols     = outcome,
    names_from  = term,
    values_from = c(edf, `ref.df`, statistic, p.value),
    names_glue  = "{term}_{.value}"
  )

# metrics
metrics_b <- sum_base_unique %>% 
  distinct(outcome, r.sq, dev.expl, n.obs, AIC, edf_total, max_concurv)

# join together
wide_base <- 
  metrics_b %>% 
  left_join(param_b,  by = "outcome") %>% 
  left_join(smooth_b, by = "outcome")

write_csv(wide_base, "Results/Gamfs/model_wide_base.csv")

# 6) SAVE BATCH PDFs ------------------------------------------------------
save_grid_pdf <- function(plot_list, filename) {
  grobs <- plot_list[sort(names(plot_list))]
  pdf(filename, width = 8.27, height = 11.69)
  grid.arrange(grobs = grobs, ncol = 5, nrow = 4)
  dev.off()
}

save_grid_pdf(dose_z_plots,    "Results/Gamfs/All_dose_z.pdf")
save_grid_pdf(dose_orig_plots, "Results/Gamfs/All_dose_orig.pdf")
save_grid_pdf(dose_zoom_plots, "Results/Gamfs/All_dose_zoom.pdf")
save_grid_pdf(int1_plots,      "Results/Gamfs/All_int1_mod1.pdf")
save_grid_pdf(int2_plots,      "Results/Gamfs/All_int1_mod2.pdf")
save_grid_pdf(int3_plots,      "Results/Gamfs/All_int1_mod3.pdf")
save_grid_pdf(
  lapply(int1_plots, function(p) p + theme(legend.position="none")),
  "Results/Gamfs/All_int1_mod_noleg.pdf"
)
save_grid_pdf(
  lapply(int2_plots, function(p) p + theme(legend.position="none")),
  "Results/Gamfs/All_int2_mod_noleg.pdf"
)
save_grid_pdf(
  lapply(int3_plots, function(p) p + theme(legend.position="none")),
  "Results/Gamfs/All_int3_mod_noleg.pdf"
)

# 7) MAP PDFs FOR SPATIAL MODELS --------------------------------------------

newdata_sub <- HEX_df_sub %>% 
  # keep only the ID
  select(H_ID) %>% 
  # bring in lon/lat from centroids_df
  left_join(centroids_df, by = "H_ID") %>% 
  # set the one‐SD treatment contrast
  mutate(N1_8Area_log_z = 1) %>% 
  # bind on covariate means for each row
  bind_cols(cov_means[rep(1, nrow(HEX_df_sub)), ])

# — 4) Find your outcomes list & create output folder —
outcomes <- c("avgAF_z","avgCHD_z","avgHF_z","avgHYP_z","avgPAD_log_z",
              "avgSTIA_z","avgOB18_z","avgCAN_z","avgDM17_z","avgNDH18_z",
              "avgCKD18_log_z","avgAST6_z","avgCOPD_z","avgRA16_z",
              "avgOST50_log_z","avgDEP_z","avgMH_log_z","avgDEM_log_z",
              "avgEP18_z","avgLD_z")
# Prepare containers for the two sets of plots
dir.create("Results/SVC_Maps/PNGs", recursive = TRUE, showWarnings = FALSE)

svc_plots_pred <- list()
svc_plots_dev  <- list()
pred_list <- list()
dev_list <- list()

for (oc in outcomes) {
  # 1) Subset to rows used in this model
  req    <- c(oc, "N1_8Area_log_z", confounders, prognostics, pca_lags)
  keep   <- complete.cases(HEX_df[, req])
  df_sub <- HEX_df[keep, ]
  W13_sub<- W13[keep, keep]
  
  # 2) Fit base GAM & compute lag‐residual
  f_base   <- as.formula(
    paste0(
      oc,
      " ~ s(N1_8Area_log_z, k=12, bs='ts') + ",
      paste(c(confounders, prognostics, pca_lags), collapse = " + ")
    )
  )
  mod_base <- bam(
    f_base,
    data     = df_sub,
    weights  = df_sub[[weights]],
    method   = "REML",
    gamma    = 1.4,
    select   = TRUE,
    discrete = TRUE
  )
  df_sub$resid_lag13 <- as.numeric(
    W13_sub %*% residuals(mod_base, type = "response")
  )
  
  # 3) Fit the SVC‐Spat model
  f_svc <- update(
    f_base,
    paste0(
      ". ~ . + resid_lag13 + ",
      "s(lon,lat,by=N1_8Area_log_z, bs='ts', k=200)"
    )
  )
  mod_svc <- bam(
    f_svc,
    data     = df_sub,
    weights  = df_sub[[weights]],
    method   = "REML",
    gamma    = 1.4,
    discrete = TRUE,
    select   = TRUE
  )
  
  saveRDS(
    mod_svc,
    file = file.path("Results/SVC_Maps", sprintf("mod_svc_%s.rds", oc))
  )
  
  # 4) Build prediction grid matching df_sub
  newdata_sub <- df_sub %>%
    select(H_ID, lon, lat) %>%
    mutate(
      N1_8Area_log_z = 1,   # one‐SD treatment
      resid_lag13    = 0    # hold residual at 0
    ) %>%
    bind_cols(cov_means[rep(1, nrow(df_sub)), ])
  
  # 5) Extract the two smooth terms
  terms_mat <- predict(mod_svc, newdata = newdata_sub, type = "terms")
  tm        <- colnames(terms_mat)
  
  treat_term <- grep("^s\\(N1_8Area_log_z\\)$", tm, value = TRUE)
  spat_term  <- grep("^s\\(lon\\s*,\\s*lat\\):N1_8Area_log_z$", tm, value = TRUE)
  
  beta_global <- terms_mat[, treat_term]
  beta_dev    <- terms_mat[, spat_term]
  
  # 6) Prepare spatial data for plotting (only those H_IDs)
  vals_pred <- tibble(
    H_ID = newdata_sub$H_ID,
    β    = beta_global + beta_dev
  )
  vals_dev  <- tibble(
    H_ID = newdata_sub$H_ID,
    β    = beta_dev
  )
  
  pred_list[[oc]] <- vals_pred %>%
    rename(beta_pred = β) %>%
    mutate(outcome = oc)
  
  dev_list[[oc]] <- vals_dev %>%
    rename(beta_dev = β) %>%
    mutate(outcome = oc)
  
  hex_pred <- hex_sp %>%
    filter(H_ID %in% vals_pred$H_ID) %>%
    left_join(vals_pred, by = "H_ID")
  
  hex_dev <- hex_sp %>%
    filter(H_ID %in% vals_dev$H_ID) %>%
    left_join(vals_dev, by = "H_ID")
  
  # 7) Compute symmetric color limits
  maxabs_pred <- max(abs(vals_pred$β), na.rm = TRUE)
  maxabs_dev  <- max(abs(vals_dev$β),  na.rm = TRUE)
  
  # 8) Build & save individual PNGs
  berry_cols <- c("#B25D91FF", "#CB87B4FF", "#EFC7E6FF",
                  "#1BB6AFFF", "#088BBEFF", "#172869FF")
  
  p_pred <- ggplot(hex_pred) +
    geom_sf(aes(fill = β), color = NA) +
    scale_fill_gradientn(
      colors   = berry_cols,
      limits    = c(-maxabs_pred, maxabs_pred),
      na.value  = "grey90",
      name      = expression(beta~"(z)")
    ) +
    ggtitle(oc) +
    theme_void() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
  
  p_dev <- ggplot(hex_dev) +
    geom_sf(aes(fill = β), color = NA) +
    scale_fill_distiller(
      type      = "div",
      palette   = "RdBu",
      direction = -1,
      limits    = c(-maxabs_dev, maxabs_dev),
      na.value  = "grey90",
      name      = expression(Delta~beta)
    ) +
    ggtitle(oc) +
    theme_void() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
  
  ggsave(
    file.path("Results/SVC_Maps/PNGs", paste0(oc, "_pred_berry.png")),
    plot   = p_pred,
    width  = 6, height = 5, dpi = 200
  )
  ggsave(
    file.path("Results/SVC_Maps/PNGs", paste0(oc, "_dev.png")),
    plot   = p_dev,
    width  = 6, height = 5, dpi = 200
  )
  
  # 9) Store for the multi‐panel PDFs
  svc_plots_pred[[oc]] <- p_pred
  svc_plots_dev[[oc]]  <- p_dev
}

all_pred <- bind_rows(pred_list)
all_dev  <- bind_rows(dev_list)

# If you want one long table with both:
all_svc <- all_pred %>%
  left_join(all_dev, by = c("H_ID","outcome"))

# Write out CSV (or RDS)
write_csv(all_svc, "Results/SVC_Maps/all_svc_predictions.csv")
# or if you prefer R binary:
saveRDS(all_svc, "Results/SVC_Maps/all_svc_predictions.rds")

# 10) Finally, write the two 20-panel PDFs

# 10a) Predicted β maps
pdf("Results/SVC_Maps/All_SVC_Predicted.pdf", width = 8.27, height = 11.69)
grid.arrange(
  grobs = svc_plots_pred,
  ncol  = 5, nrow = 4,
  top   = textGrob(
    "Predicted standardized spatial coefficients β(lon,lat)",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()

# 10b) Δβ maps
pdf("Results/SVC_Maps/All_SVC_Deviation.pdf", width = 8.27, height = 11.69)
grid.arrange(
  grobs = svc_plots_dev,
  ncol  = 5, nrow = 4,
  top   = textGrob(
    "Standardized spatial deviation Δβ(lon,lat)",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()
