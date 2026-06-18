########################################################################
#### Externalizing — Prepare Data for Mplus format for Twin Model   #### 
#### 01.04.2026 - 11.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

library(dplyr)
library(tidyr)

# ============================================================
# 1. DATA LOADING 
# ============================================================

# Define path for data files
path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data")

# Load preprocessed data
load(file.path(path, "04_twinlife_externalizing_cfa.rda"))

head(df_fs)

# ============================================================
# 2. DATA FORMAT PREP FOR MPLUS
# ============================================================

# Prepare df to export PC1 in wide format with row per fid for Mplus
# Also export EXT and first order factor scores for heritability sensitivity analysis
df_pc1_wide <- df_fs %>%
  group_by(fid) %>%
  filter(n() == 2) %>% # complete twin pairs only, no single twin participants
  arrange(pid, .by_group = TRUE) %>%
  mutate(twin = paste0("t", row_number())) %>%
  ungroup() %>%
  select(fid, pid, sex, zygosity, age_yrs, age_group, PC1, EXT, hyp, att, srg, con, twin) %>%
  pivot_wider(
    id_cols = c(fid, zygosity, age_yrs, age_group),
    names_from = twin,
    values_from = c(pid, sex, PC1, EXT, hyp, att, srg, con),
    names_sep = "_"
  ) %>%
  select(
    fid,
    pid_t1, pid_t2,
    sex_t1, sex_t2,
    zygosity,
    age_yrs,
    age_group,
    PC1_t1, PC1_t2,
    EXT_t1, EXT_t2,
    hyp_t1, hyp_t2,
    att_t1, att_t2,
    srg_t1, srg_t2,
    con_t1, con_t2
  )

head(df_pc1_wide)


# ============================================================
# 3. AGE CONSISTENCY
# ============================================================

# Check number of twin pairs per age
nrow(df_pc1_wide) # more twin pair rows than in initial dataset ?

# Count fids
df_pc1_wide %>%
  count(fid) %>%
  count(n) # correct: separate rows for 11 twins with no common age

# Exclude age inconsistent twin pairs for twin model
df_pc1_wide <- df_pc1_wide %>%
  group_by(fid) %>%
  filter(n() == 1) %>%
  ungroup()

nrow(df_pc1_wide) # final sample size: 3685 twin pairs

# Check for opposite gender twins to great gender-zygosity groups for Mplus
df_pc1_wide %>%
  filter(sex_t1 != sex_t2) %>%
  select(fid, pid_t1, pid_t2, sex_t1, sex_t2) # none


# ============================================================
# 4. CREATE AGE, GENDER, ZYGOSITY GROUPS FOR MPLUS
# ============================================================

# Create age group- and gender-zygosity groups for Mplus
df_pc1_wide_mplus <- df_pc1_wide %>%
  mutate(
    # numeric age group
    age_group = as.integer(age_group),
    # zygosity-gender group variable
    zyg_gender = case_when(
      zygosity == 1 & sex_t1 == 1 ~ 1, # 1 = MZ male
      zygosity == 2 & sex_t1 == 1 ~ 2, # 2 = DZ male
      zygosity == 1 & sex_t1 == 2 ~ 3, # 3 = MZ female
      zygosity == 2 & sex_t1 == 2 ~ 4, #4 = DZ female
      TRUE ~ NA_real_
    ),
    
    # zygosity-age-group variable
    zyg_age = case_when(
      age_group == 1 & zygosity == 1 ~ 1,  #   1  = MZ age group 1
      age_group == 1 & zygosity == 2 ~ 2,  #   2  = DZ age group 1
      age_group == 2 & zygosity == 1 ~ 3,  #   3  = MZ age group 2
      age_group == 2 & zygosity == 2 ~ 4,  #   4  = DZ age group 2
      age_group == 3 & zygosity == 1 ~ 5,  #   5  = MZ age group 3
      age_group == 3 & zygosity == 2 ~ 6,  #   6  = DZ age group 3
      age_group == 4 & zygosity == 1 ~ 7,  #   7  = MZ age group 4
      age_group == 4 & zygosity == 2 ~ 8,  #   8  = DZ age group 4
      age_group == 5 & zygosity == 1 ~ 9,  #   9  = MZ age group 5
      age_group == 5 & zygosity == 2 ~ 10, #   10 = DZ age group 5
      ),
      
      # age-group-gender-zygosity variable
      age_gender_zyg = case_when(
        age_group == 1 & sex_t1 == 1 & zygosity == 1 ~ 1,  # Male MZ age group 1 = age_gender_zyg group 1 etc.
        age_group == 1 & sex_t1 == 1 & zygosity == 2 ~ 2,  # Male DZ age group 1 = age_gender_zyg group 1 etc.
        age_group == 1 & sex_t1 == 2 & zygosity == 1 ~ 3,  # Female MZ age group 1 = age_gender_zyg group 3 etc.
        age_group == 1 & sex_t1 == 2 & zygosity == 2 ~ 4,  # Female DZ age group 1 = age_gender_zyg group 3 etc.
        
        age_group == 2 & sex_t1 == 1 & zygosity == 1 ~ 5,
        age_group == 2 & sex_t1 == 1 & zygosity == 2 ~ 6,
        age_group == 2 & sex_t1 == 2 & zygosity == 1 ~ 7,
        age_group == 2 & sex_t1 == 2 & zygosity == 2 ~ 8,
        
        age_group == 3 & sex_t1 == 1 & zygosity == 1 ~ 9,
        age_group == 3 & sex_t1 == 1 & zygosity == 2 ~ 10,
        age_group == 3 & sex_t1 == 2 & zygosity == 1 ~ 11,
        age_group == 3 & sex_t1 == 2 & zygosity == 2 ~ 12,
        
        age_group == 4 & sex_t1 == 1 & zygosity == 1 ~ 13,
        age_group == 4 & sex_t1 == 1 & zygosity == 2 ~ 14,
        age_group == 4 & sex_t1 == 2 & zygosity == 1 ~ 15,
        age_group == 4 & sex_t1 == 2 & zygosity == 2 ~ 16,
        
        age_group == 5 & sex_t1 == 1 & zygosity == 1 ~ 17,
        age_group == 5 & sex_t1 == 1 & zygosity == 2 ~ 18,
        age_group == 5 & sex_t1 == 2 & zygosity == 1 ~ 19,
        age_group == 5 & sex_t1 == 2 & zygosity == 2 ~ 20,
        
        TRUE ~ NA_real_
        )
  )

# check if all rows were assigned
df_pc1_wide_mplus %>%
  filter(is.na(zyg_age) | is.na(zyg_gender)) # yes

head(df_pc1_wide_mplus)

# Remove twins without PC1 to retain final sample
df_pc1_wide_mplus <- df_pc1_wide_mplus %>%
  filter(!is.na(PC1_t1), !is.na(PC1_t2))

# Final sample size number of twin pairs with PC1
length(unique(df_pc1_wide_mplus$fid)) #3389


# ============================================================
# 5. MZ DZ TWIN CORRELATIONS
# ============================================================

# Remove twins without PC1 to retain final sample
df_pc1_wide_mplus <- df_pc1_wide_mplus %>%
  filter(!is.na(PC1_t1), !is.na(PC1_t2))

pc1_twin_corr <- df_pc1_wide_mplus %>%
  group_by(zygosity) %>%
  summarise(
    n = n(),
    twin_correlation = cor(PC1_t1, PC1_t2)
  )

pc1_twin_corr

factor_twin_corr <- df_pc1_wide_mplus %>%
  group_by(zygosity) %>%
  summarise(
    n_hyp = sum(complete.cases(hyp_t1, hyp_t2)),
    r_hyp = cor(hyp_t1, hyp_t2, use = "complete.obs"),
    
    n_att = sum(complete.cases(att_t1, att_t2)),
    r_att = cor(att_t1, att_t2, use = "complete.obs"),
    
    n_srg = sum(complete.cases(srg_t1, srg_t2)),
    r_srg = cor(srg_t1, srg_t2, use = "complete.obs"),
    
    n_con = sum(complete.cases(con_t1, con_t2)),
    r_con = cor(con_t1, con_t2, use = "complete.obs")
  )

factor_twin_corr

# ============================================================
# 6. SAMPLE SIZE PER AGE GROUP
# ============================================================

sample_sizes <- df_pc1_wide_mplus %>%
  count(age_group, name = "n_pairs") %>%
  arrange(age_group)

sample_sizes


# Per zygosity
sample_sizes_zyg <- df_pc1_wide_mplus %>%
  count(age_group, zygosity, name = "n_pairs") %>%
  arrange(age_group, zygosity)

sample_sizes_zyg

# ============================================================
# 7. EXPORT DATA
# ============================================================

# Export for Mplus Twin Model
write.table(
  df_pc1_wide_mplus,
  file = file.path(path, "05_twinlife_externalizing_formplus.dat"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t",
  na = "." # create missingness marker for Mplus
)


