########################################################################
#### Externalizing — Factor Invariance Testing in TwinLife          #### 
#### 01.04.2026 - 11.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

library(lavaan)

# ============================================================
# 1. DATA LOADING 
# ============================================================

# Define path for data files
path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data")

# Define path for factor model
path_model <- ("/Volumes/MPRG-Biosocial/Projects/04_data_analysis/012_SHIP_BASE_TwinLife_SOEP/Externalizing_Genomics/models")

# Load preprocessed data
load(file.path(path, "04_twinlife_externalizing_model_data.rda"))

# ============================================================
# 2. INVARIANCE TESTING ACROSS AGE GROUPS
# ============================================================

# Age in years as grouping variable
df_model$age_group <- as.factor(df_model$age_yrs_pq)

# Test invariance across age groups (can't converge with age continously)
# Define groups (tried groups of 2 yrs, didn't converge -> 3 yrs worked!)
df_model$age_group <- cut(df_model$age_yrs_pq, 
                       breaks = c(10, 15, 21, 26), 
                       labels = c("11-15",  # 11-15 early to middle adolescence 
                                  "16-21",  # 16-21 late adolescence to early emerging adulthood 
                                  "22-26")) # 22-26 later emerging adulthood

# Load final model
source(file.path(path_model, "model_final.R"))

# 1. Fit configural model
fit_configural <- cfa(
  model = model_final, 
  data = df_model, 
  group = "age_group", # fit factor separately across groups
  estimator = "MLR", 
  missing = "FIML",
  cluster = "fid",
  
  # Use the first indicator to set the scale for stability
  std.lv = FALSE
)

summary(fit_configural, fit.measures = TRUE)


# 2. Test metric invariance 
fit_metric <- cfa(
  model = model_final, 
  data = df_model, 
  group = "age_group",
  group.equal = "loadings", # fix loadings to be equal across age groups
  estimator = "MLR", 
  missing = "FIML", 
  cluster = "fid"
)

summary(fit_metric, fit.measures = TRUE)

# Define metrics to compare for both models models
invariance_metrics = c("cfi", "cfi.scaled", "cfi.robust", "rmsea", "rmsea.scaled", "rmsea.robust", "srmr")

# Create table for invariance test results
invariance_results <- rbind(
  Configural = fitMeasures(
    fit_configural,
    invariance_metrics
  ),
  Metric = fitMeasures(
    fit_metric,
    invariance_metrics
  )
)

invariance_results <- round(invariance_results, 3)

# add difference row 
invariance_results <- rbind(
  invariance_results,
  Delta = round(
    invariance_results["Metric", ] - invariance_results["Configural", ],
    3
  )
)

print(invariance_results) # metrics within acceptable ranges for both configural and metric invariance


# ============================================================
# 2. INVARIANCE TESTING ACROSS GENDER
# ============================================================

# 1. Configural invariance
fit_configural_sex <- cfa(
  model = model_final,
  data = df_model,
  group = "sex",
  estimator = "MLR",
  missing = "FIML",
  cluster = "fid",
  std.lv = FALSE
)

summary(
  fit_configural_sex,
  fit.measures = TRUE,
  standardized = TRUE
)

# 2. Metric invariance
fit_metric_sex <- cfa(
  model = model_final,
  data = df_model,
  group = "sex",
  group.equal = "loadings",
  estimator = "MLR",
  missing = "FIML",
  cluster = "fid",
  std.lv = FALSE
)

summary(
  fit_metric_sex,
  fit.measures = TRUE,
  standardized = TRUE
)


# Create table for invariance test results
invariance_results_sex <- rbind(
  Configural = fitMeasures(
    fit_configural_sex,
    invariance_metrics
  ),
  Metric = fitMeasures(
    fit_metric_sex,
    invariance_metrics
  )
)

# Add difference row 
invariance_results_sex <- rbind(
  invariance_results_sex,
  Delta =
    invariance_results_sex["Metric", ] -
    invariance_results_sex["Configural", ]
)

invariance_results_sex <- round(
  invariance_results_sex,
  3
)

print(invariance_results_sex)
