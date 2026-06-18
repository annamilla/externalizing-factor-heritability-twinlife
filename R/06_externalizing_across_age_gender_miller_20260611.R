########################################################################
#### Externalizing — EXT variation across age and gender            #### 
#### 14.04.2026 - 14.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

library(dplyr)   
library(lmtest)
library(sandwich)

# ============================================================
# 1. DATA LOADING 
# ============================================================

path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data")

# Load final model data
load(file.path(path, "04_twinlife_externalizing_cfa.rda"))

head(df_fs)


# ============================================================
# 2. ADD AGE GROUPS
# ============================================================


# Define age groups according to invariance testing
df_fs$age_group <- cut(df_fs$age_yrs, 
                       breaks = c(10, 15, 21, 26), 
                       labels = c("11-15",  # 11-15 early to middle adolescence 
                                  "16-21",  # 16-21 late adolescence to early emerging adulthood 
                                  "22-26")) # 22-26 later emerging adulthood

head(df_fs)

# age group as numeric
df_fs <- df_fs %>%
  mutate(
    age_group = as.numeric(factor(
      age_group,
      levels = c("11-13", "14-16", "17-19", "20-22", "23-26"))
    ))


# ============================================================
# 3. REGRESSION PREPARATION GENDER AND AGE CODING
# ============================================================

# Create gender effect-coded, age centered and interaction term in df
df_fs <- df_fs %>%
  mutate(
    gender = ifelse(sex == 2, 0.5, -0.5), # effect code sex female as +0.5, male as -0.5
    age_group_centered = age_group - mean(age_group, na.rm = TRUE),
    age_x_gender = age_group_centered * gender
  )

head(df_fs)


# ============================================================
# 4. REGRESS PC1 ON AGE, GENDER, INTERACTION
# ============================================================

# 1. Regression model with PC1 as the outcome and age_group_centered, gender, and their interaction as predictors
lm_pc1 <- lm(
  PC1 ~ age_group_centered * gender,
  data = df_fs
)

# Coefficient table clustered by twin pairs (fid)
pc1_results <- lmtest::coeftest(
  lm_pc1,
  vcov = sandwich::vcovCL(
    lm_pc1,
    cluster = df_fs$fid
  )
)

# Convert to df and apply FDR correction
pc1_table <- data.frame(
  term = rownames(pc1_results),
  estimate = pc1_results[, 1],
  se = pc1_results[, 2],
  t_value = pc1_results[, 3],
  p_value = pc1_results[, 4],
  row.names = NULL
) %>%
  mutate(
    p_fdr = p_value,
    p_fdr = replace(
      p_fdr,
      term != "(Intercept)", # exclude intercept from multiple testing
      p.adjust(p_value[term != "(Intercept)"], method = "fdr")
    )
  )

pc1_table   


# 2. Test for non-linear age effects
# Make age group categorical
df_fs <- df_fs %>%
  mutate(
    age_group_f = factor(
      age_group,
      levels = 1:3,
      labels = c("11-15", "16-21", "22-26")
    )
  )

# Non-linear categorical age model
lm_pc1_age_cat <- lm(
  PC1 ~ age_group_f * gender,
  data = df_fs
)

# Coefficients table
pc1_age_cat_results <- lmtest::coeftest(
  lm_pc1_age_cat,
  vcov = sandwich::vcovCL(lm_pc1_age_cat, cluster = df_fs$fid)
)

# Convert to df and FDR correct
pc1_age_cat_table <- data.frame(
  term = rownames(pc1_age_cat_results),
  estimate = pc1_age_cat_results[, 1],
  se = pc1_age_cat_results[, 2],
  t_value = pc1_age_cat_results[, 3],
  p_value = pc1_age_cat_results[, 4],
  row.names = NULL
) %>%
  # add column for FDR corrected p value, 
  mutate(
    p_fdr = p_value,
    p_fdr = replace(
      p_fdr,
      term != "(Intercept)", # exclude intercept from multiple testing
      p.adjust(p_value[term != "(Intercept)"], method = "fdr")
    )
  )

pc1_age_cat_table # roughly linear!


# ============================================================
# 5. REGRESS FIRST ORDER FACTORS ON AGE, GENDER, INTERACTION
# ============================================================

# Define first order factor columns
fo_fs <- c("hyp", "att", "srg", "con")

# Test for age, gender variation across first order factor scores
lm_results <- lapply(fo_fs, function(outcome) {
  model <- lm(
    as.formula(paste0(outcome, " ~ age_group_centered * gender")),
    data = df_fs
  )
  # coefficients 
  lmtest::coeftest(
    model,
    vcov = sandwich::vcovCL(model, cluster = df_fs$fid)
  )
})

# Add first order factor names to result table
names(lm_results) <- fo_fs

# Create regression table
lm_table <- do.call(
  rbind,
  lapply(names(lm_results), function(outcome) {
    
    x <- lm_results[[outcome]]
    
    data.frame(
      outcome = outcome,
      term = rownames(x),
      estimate = x[, 1],
      se = x[, 2],
      t_value = x[, 3],
      p_value = x[, 4],
      row.names = NULL
    )
  })
)

lm_table

# FDR correction  
lm_table_fdr <- lm_table %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  # add column for fdr corrected p value
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr")
  ) %>%
  ungroup()


# Table with significant effects only
fo_results_sig <- lm_table_fdr %>%
  filter(p_value < .05)

fo_results_sig


# ============================================================
# 6. EXPORT DATA
# ============================================================

# Export data including age group column overwriting previously created file
write.csv(df_fs, file.path(path, "04_twinlife_externalizing_cfa.csv"), row.names = FALSE)
save(df_fs, file = file.path(path, "04_twinlife_externalizing_cfa.rda"))

