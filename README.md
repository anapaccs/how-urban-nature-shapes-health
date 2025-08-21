# how-urban-nature-shapes-health
This repository contains the datasets, results, and scripts that accompany the PhD thesis:

**"How Urban Nature Shapes Health: Unravelling the Spatial Complexity of Nature–Health Systems"**  

Department of Land Economy

University of Cambridge

August 2025  

Author: Ana Paula Albuquerque Campos Costalonga Seraphim  

Each file or folder is linked to the thesis, with references to the relevant **chapter, section, and page numbers** where it is introduced or discussed.  

---

## 📊 Data Tables

- **PolicyReports-DHs-HealthOutcomes-PriorityThemes.xlsx**  
  Matrix indicating which policy reports, reviewed in the Chapter 2 scoping study, reference each Determinant of Health (DH), Health Outcome, and Priority Theme.  

  *Referenced in Chapter 2, Section 2.4, p. 21.*

- **Bibliometric-Selected-Articles-Natural-Spaces-Health.xlsx**  
  Dataset containing bibliographic information of the articles included in Chapter 3's bibliometric review. Fields include: RefID, Author, Title, Abstract, Keywords, DOI, ISSN, Journal, Issue, Volume, Pages, Place Published, Publisher, Type, and Year.

  *Referenced in Chapter 3, Section 3.2.2, p. 60.*

- **Cooccurrence-Network-Articles-Title-Abstract-Terms.xlsx**  
  Dataset containing the terms resulting from the co-occurrence network analysis of selected articles' titles and abstracts in the bibliometric review (Chapter 3). Each term is provided with associated details, including: ID, Category, Period, x, y, Cluster, Weight (Links), Weight (Total Link Strength), Weight (Occurrences), and Score (Average Publication Year). Terms are further organised into time intervals.

  *Referenced in Chapter 3, Section 3.2.2, p. 60.*

- **Morans-I-Results.xlsx**  
  Results of incremental Moran’s I across 30 distance bands (from 350 m up to 12 km) using a Euclidean approach for indicators of Chapter 4's spatial dataset. Fields reported: var, band, lower, upper, midpoint, I, Z, p_value, effect.

  *Referenced in Chapter 4, Section 4.3.2.2, p. 101.*

  *Referenced in Chapter 4, Section 4.4.2, p. 112.*

  *Referenced in Chapter 5, Section 5.2.5.1, p. 129.*

- **Correlation-Matrix-Table.xlsx**  
  Pairwise Spearman correlation matrix summarising relationships among the 46 transformed variables in Chapter 5 Analysis, providing an initial exploration of broad co-variation patterns. 

  *Referenced in Chapter 5, Section 5.2.4.1, p. 125.*

---

## 🗺️ Spatial Data Framework (`/spatial_data_framework`)
Complete spatial dataset developed for the Greater London case study in Chapter 4, provided in both CSV and GeoPackage formats, together with metadata.

 *Referenced in Chapter 4, Section 4.4.1, p. 108.*
 
- **London_HEX_NH.csv**  
  Hexagon-level dataset (350m resolution) containing variables on nature, socioeconomic, environmental, and health indicators.

- **London_HEX_NH.gpkg**  
  GeoPackage file with the same dataset for spatial analysis in GIS.

- **Metadata_London_HEX_NH.xlsx**  
  Metadata describing all variables in the spatial dataset.

- **Summary-Statistics-Indicators.xlsx**  
 Descriptive statistics of all indicators included in the spatial dataset. It includes: indicator, n, nmiss, min, Q1, median, mean, Q3, max, sd, cv, skew, kurt, IQR, mad, range, and outlier.

*Referenced in Chapter 4, Section 4.4.1, p. 108.*

*Referenced in Chapter 5, Section 5.2.3, p. 124.*

---

## 💻 R Scripts (`/analysis_R_scripts`)
This folder contains the analytical workflow implemented in R (v4.2.0) to quantify the relationship between natural space exposure and health outcomes. The scripts cover spatial lag computation, multicollinearity testing, Generalised Propensity Score estimation with stabilised IPW, and spatial Generalised Additive Models (GAMs) for non-linear dose–response analysis, moderation, and spatial heterogeneity.  

*Referenced in Chapter 5, Section 5.2.5, p. 129.*

- **A1.Lags_PCA.R**  
  Script for generating spatial lags and running Principal Component Analysis.  

- **A2.VIF.R**  
  Script for Variance Inflation Factor analysis to test multicollinearity.  

- **A3.GPS_cof.R**  
  Script for estimating Generalised Propensity Scores with confounders.
  
- **A3.GPS_cov.R**  
  Script for estimating Generalised Propensity Scores with covariates.  
 
- **A4.Spatial_GAM.R**  
  Script for running spatial Generalised Additive Models with GPS weighting.  

- **A5.Lagrange.R**  
  Script for Lagrange Multiplier tests for spatial dependence.  
 
---

## 🌀 Spatial GAM + GPS Results (`/spatial_GPS-GAM_results`)
This folder contains the results of our Spatial Generalised Additive Model (GAM) with Generalised Propensity Score (GPS) weighting models applied to twenty health outcomes across Greater London. The outputs include model summaries (Excel files) and full-size plots illustrating non-linear dose–response relationships and ecological moderation.  
 
 *Referenced in Chapter 5, Section 5.3, p. 137.*
 
- **Base-Model-Results.xlsx**  
  Summary results of the baseline GPS–GAM model.  
 
- **Covariate-lag-Model-Results.xlsx**  
  Summary results of GPS–GAM models, including covariate spatial lags.   

- **Residual-lag-Model-Results.xlsx**  
  Summary results of GPS–GAM models, including both covariate and residual spatial lags (final specification).  

### Subfolders by Model Type
- **covariate-lag_plots/**  
  Full plots for the covariate-lag GPS–GAM models (20 health outcomes). Each outcome includes:   
  - Dose–response on the z-scale  
  - Dose–response on the original scale  
  - Dose–response on the original scale zoomed for Nature Expose dose between 0 and 250 hectares 
  - GAM-check residual spatial plots
  
- **residual-lag_plots/**  
  Full plots for the residual-lag GPS–GAM models (20 health outcomes). Each outcome includes:  
  - Dose–response on the z-scale  
  - Dose–response on the original scale  
  - Dose–response on the original scale zoomed for Nature Expose dose between 0 and 250 hectares
  - GAM-check residual spatial plots  
 
- **residual-lag_moderation_plots/**  
  Moderation plots for the residual-lag GPS–GAM models (20 health outcomes). Each outcome includes dose–response curves on the z-scale, moderated by:   
  - Biodiversity index
  - Tree cover  
  - Vegetation cover   

## 📌 Notes 
- All results are reproducible using the R scripts provided (`/analysis_R_scripts`), together with the spatial data framework (`/spatial_data_framework`).  

---

