########################################################################
#### Externalizing — Phenotypic data cleaning in TwinLife           #### 
#### 25.02.2026 - 19.05.2026                                        ####
#### Clean and prepare TwinLife externalizing data                  ####
#### main Author: Anna Miller                                       #### 
########################################################################

# ============================================================
# 0. IMPORT LIBRARIES FOR DATA MANIPULATION AND PLOT
# ============================================================

library(tidyr)
library(haven)
library(dplyr)
library(ggplot2)


# ============================================================
# 1. DATA LOADING
# ============================================================

# Define path
path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data")

# Load twins self-reported externalizing dataset
df_twinlife <- load(file.path(path, "01_twinLife_externalizing_data.rda"))


# ============================================================
# 2. MISSING DATA
# ============================================================


# Get all unique values across df_ext to check missing codes with TwinLife data documentation
ext_values <- df_ext %>%
  select(-pid) %>%                      # exclude pid column
  mutate(across(everything(), as.numeric)) %>%  
  unlist() %>%                          
  unique() %>%                          
  sort()                                

head(ext_values, 20) 

# Define missing codes from previous script and data documentation
missing_codes <- c(-99, -98, -95, -94, -93, -92, -90, -89, -87)

# Missing code 81 (haven't had sexual intercourse yet):
# provides valuable data on whether someone ever had sex at that age!
# Create variable for ever had sex
df_ext <- df_ext %>%
  mutate(
    sex_ever = if_else(!sex_prot %in% missing_codes, 
                       if_else(sex_prot == -81, 0, 1), 
                       NA))
# Check variable values
table(df_ext$sex_ever)

# Redefine data missing codes for externalizing variables
missing_codes <- c(-99, -98, -95, -94, -93, -92, -90, -87, -86, -85, -81)

# Replace missing_codes values in all variables with NA
df_ext <- df_ext %>%
  mutate(across(everything(), ~ ifelse(. %in% missing_codes, NA, .)))

# ============================================================
# 3. DUPLICATE VARIABLES
# ============================================================

# Inspect and harmomonize duplicate variables
# 1. Smok and smok_ever: Look at rows with available data
df_smok <- df_ext %>%
  filter(!is.na(smok) |!is.na(smok_ever)) %>%
  select(pid, smok, smok_ever) # Looks like either smok_init or smok_init2 are available across rows

# Check whether both smok and smok_ever are available in any row before merging
df_ext %>%
  filter(!is.na(smok) & !is.na(smok_ever)) %>%
  select(pid, smok, smok_ever) # no, exclusive

# Harmonize smok and smok_ever
df_ext <- df_ext %>%
  mutate(
    # Recode smok_ever
    smok_ever_recoded = case_when(
      smok_ever == 1 ~ 2, # 1 (yes) -> 2 (smoker)
      smok_ever == 2 ~ 6, # 2 (no) -> 6 (no)
      TRUE ~ NA_real_
    ),
    # Merge prioritizing smok (more detailed), then recoded smok_ever
    smok = coalesce(as.numeric(smok), smok_ever_recoded)
  ) %>%
  select(-smok_ever_recoded, -smok_ever) # drop the redundant variables


# 2. Smoking initiation: look at rows with available data
df_smok_init <- df_ext %>%
  filter(!is.na(smok_init) |!is.na(smok_init2)) %>%
  select(pid, smok_init, smok_init2) 

# Confirm whether both smok_init and smok_init2 are available in any row before merge
df_ext %>%
  filter(!is.na(smok_init) & !is.na(smok_init2)) %>%
  select(pid, smok_init, smok_init2) # no, exclusive

# Merge smok_init variables
df_ext <- df_ext %>%
  mutate(smok_init3 = ifelse(!is.na(smok_init), smok_init, smok_init2)) # create third variable to check

# Check merge
df_ext %>%
  summarise(
    smok_init = sum(!is.na(smok_init)),
    smok_init2 = sum(!is.na(smok_init2)),
    smok_init3 = sum(!is.na(smok_init3)),
    diff = smok_init3 - smok_init2 - smok_init # completeness check
  )

# Keep harmonized smok_init variable only
df_ext <- df_ext %>%
  mutate(smok_init = smok_init3) %>% # overwrite smok_init with harmonized variable
  select(-smok_init2, -smok_init3) # drop other smok_init variables


# 4. Smoking alone
# Confirm whether both smok_solo and smok_solo2 are available in any row
df_ext %>%
  filter(!is.na(smok_solo) & !is.na(smok_solo2)) %>%
  select(pid, smok_solo, smok_solo2) # no, exclusive

# Merge smok_alone variables
df_ext <- df_ext %>%
  mutate(smok_solo = ifelse(!is.na(smok_solo), smok_solo, smok_solo2)) %>%
  select(-smok_solo2)  # drop other smok_solo variable


# 5. Alcohol initiation
# Confirm whether both alc_init and alc_init2 are available in any row
df_ext %>%
  filter(!is.na(alc_init) & !is.na(alc_init2)) %>%
  select(pid, alc_init, alc_init2) # no, exclusive

# 6. Merge alc_init variables
df_ext <- df_ext %>%
  mutate(alc_init = ifelse(!is.na(alc_init), alc_init, alc_init2)) %>%
  select(-alc_init2) # Drop alc_init2 


# 7. Alcohol frequency excessive
# Confirm whether both alc_freq_exc and alc_freq_exc2 are available in any row
df_ext %>%
  filter(!is.na(alc_freq_exc) & !is.na(alc_freq_exc2)) %>%
  select(pid, alc_freq_exc, alc_freq_exc2) # no, exclusive

#  Merge alc_freq_exc variables
df_ext <- df_ext %>%
  mutate(alc_freq_exc = ifelse(!is.na(alc_freq_exc), alc_freq_exc, alc_freq_exc2)) %>%
  select(-alc_freq_exc2) # Drop alc_freq_exc2


# 8. devi_anger: look at rows with available data
df_devi_anger <- df_ext %>%
  filter(!is.na(devi_anger)  | !is.na(devi_anger2)) %>%
  select(pid, devi_anger, devi_anger2) 


# Confirm whether both devi_anger and devi_anger 2 are available in any row
df_ext %>%
  filter(!is.na(devi_anger) & !is.na(devi_anger2)) %>%
  select(pid, devi_anger, devi_anger2) # no, exclusive

# Merge devi_anger  variables
df_ext <- df_ext %>%
  mutate(devi_anger = ifelse(!is.na(devi_anger), devi_anger, devi_anger2)) %>%
  select(-devi_anger2)  # Drop devi_anger2 


# 9. alc_n_beer_exc: look at rows with available data
df_alc_n_beer_exc <- df_ext %>%
  filter(!is.na(alc_n_beer_exc) | !is.na(alc_n_beer_exc2)) %>%
  select(pid, alc_n_beer_exc, alc_n_beer_exc2) 

# Confirm whether both alc_n_beer_exc and alc_n_beer_exc2 are available in any row
df_ext %>%
  filter(!is.na(alc_n_beer_exc) & !is.na(alc_n_beer_exc2)) %>%
  select(pid, alc_n_beer_exc, alc_n_beer_exc2) # no, exclusive

# 10. Merge alc_n_beer_exc  variables
df_ext <- df_ext %>%
  mutate(alc_n_beer_exc = ifelse(!is.na(alc_n_beer_exc), alc_n_beer_exc, alc_n_beer_exc2)) %>%
  select(-alc_n_beer_exc2) # Drop alc_n_beer_exc2


# ============================================================
# 4. CLEAN AGE AT DATA COLLECTION
# ============================================================

# Create variable for age at time of individual questionnaire (pq)
# First check if any rows where: year/month of pq is available but age/year/month of family questionnaire (fq) is missing
df_check <- df_ext %>%
  filter(
    !is.na(yea_pq) & !is.na(mon_pq) &                               # year/month of pq available
      (is.na(age_mon_fq) | is.na(age_yrs_fq) | is.na(age_mon_fq) |  # age at family questionnaire (fq) missing
         is.na(yea_fq) | is.na(mon_fq))                             # year/month of family questionnaire (fq) missing
  )

nrow(df_check) # no missing year/month of fq for all rows with year/month of pq 
                # -> can use fq age to compute pq age!

# Compute age in years and months at pq
df_ext <- df_ext %>%
  mutate(
    # compute time difference between pq and fq
    diff_mon = (yea_pq - yea_fq) * 12 + (mon_pq - mon_fq),
    diff_yrs = yea_pq - yea_fq,
    
    # compute age at pq where pq year and month are not missing (-95 according to data documentation)
    age_mon_pq = if_else(mon_pq != -95, age_mon_fq + diff_mon, NA),
    age_yrs_pq = if_else(yea_pq != -95, age_yrs_fq + diff_yrs, NA)
  ) %>%
  select(-diff_mon, -diff_yrs) # drop diff variables again

# check correctness of age at pq computations where year/month of fq don't match year/month of pq
df_check <- df_ext %>%
  filter(
    (yea_fq != yea_pq) | (mon_fq != mon_pq)
  ) %>%
  select(pid, birth_month, birth_year, yea_fq, mon_fq, yea_pq, mon_pq, age_yrs_fq, age_mon_fq, age_yrs_pq, age_mon_pq, age_twins_yrs_fq)

# Check some random samples
View(df_check)

# Check age range at time of individual questionnaire
range(df_ext$age_yrs_pq, na.rm = TRUE) # as expected according to data documentation
range(df_ext$age_mon_pq, na.rm = TRUE) # as expected according to data documentation

# Check birth year range by examining age at last questionnaire (for study description)
df_birth_year <- df_ext %>%
  select(pid, fid, age_yrs_pq, yea_pq) %>%
  filter(
    !is.na(age_yrs_pq )) %>%
  mutate(
  birth_year = yea_pq - age_yrs_pq    # compute birth year from pq year and age
  )

head(df_birth_year)

min(df_birth_year$birth_year)
max(df_birth_year$birth_year)


# Drop variables from df_ext that are not needed anymore
df_ext <- df_ext %>%
  select(-yea_fq, -mon_fq, -age_yrs_fq, -age_mon_fq, -age_twins_yrs_fq, -mon_pq, -yea_pq, -birth_month, -birth_year, -cgr, -wid)

# Sanity check all rows on one example pid
df_check <- df_ext %>% filter(pid == 110445002)
head(df_check) 
View(df_check) # multiple entries for different timepoints, as expected


# ============================================================
# 5. VARIABLE DISTRIBUTION ACROSS AGE
# ============================================================

# Count pids with available data per variable
df_ext_counts <- df_ext %>%
  mutate(across(-pid, ~ as.numeric(.x))) %>%  # convert to numeric
  pivot_longer(cols = -pid, names_to = "variable", values_to = "value") %>%  # long format
  filter(!is.na(value) & !(value %in% missing_codes)) %>%  
  group_by(variable) %>% # count per variable
  summarise(n_unique_pids = n_distinct(pid)) %>%  # count unique pids
  arrange(desc(n_unique_pids)) 

View(df_ext_counts) # completeness for age, sex, fid confirmed (4093 twins with zygosity data)

# Create df grouped by variable category
df_ext_counts_grouped <- df_ext_counts %>%
  arrange(variable)

View(df_ext_counts_grouped)


# Check coverage of all variables across age groups, 
# to select variables and/or age groups for analyses
# Define demographic variables other than age in years to deselect in age coverage check
demo_vars <- c("pid", "fid", "sex", "zygosity", "age_mon_pq")

# Calculate completeness (0 to 1) for every variable by age group
age_map <- df_ext %>%
  select(age_yrs_pq, everything(), -all_of(demo_vars)) %>%
  group_by(age_yrs_pq) %>%
  summarise(across(everything(), ~ sum(!is.na(.)) / n())) %>%
  pivot_longer(-age_yrs_pq, names_to = "variable", values_to = "completeness") %>%
  pivot_wider(names_from = age_yrs_pq, values_from = completeness)

# View the full matrix
View(age_map)

# Save
write.csv(age_map, file.path(path, "02_twinLife_externalizing_age_map.csv"), row.names = FALSE)

# Create heatmap
# Convert data to long format is for ggplot
plot_data <- df_ext %>%
  group_by(age_yrs_pq) %>%
  summarise(across(everything(), ~ mean(!is.na(.x)))) %>%
  pivot_longer(-age_yrs_pq, names_to = "variable", values_to = "completeness") %>%
  filter(!variable %in% demo_vars) %>% # filter out demographic variables 
  arrange(variable) %>% # sort by name and convert to factor for sorting
  mutate(variable = factor(variable, levels = rev(unique(variable))))

# Create heatmap
ggplot(plot_data, aes(x = factor(age_yrs_pq), y = variable, fill = completeness)) +
  # create tiles
  geom_tile(color = "white", size = 0.1) + 
  # add labels, format values to 2 decimal places
  geom_text(aes(label = sprintf("%.2f", completeness)), size = 2, color = "black") + 
  # fix cell scale
  coord_fixed() +
  # style
  scale_fill_gradientn(
    colors = c("#FFFFE0", "#99D8C9", "#081D58"),  # set color gradient
    values = c(0, 0.5, 1),
    limits = c(0, 1),
    name = "Proportion Complete"   # how many participants have observations out of total number of participants in that age group
  ) +
  # Format title and axes
  labs(
    title = "Coverage of Variables Across Age (Years)", 
    x = "Age (Years)",
    y = "Variable Name"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    panel.grid = element_blank(), # remove background grid
    legend.position = "right",
    legend.key.height = unit(2, "cm") # makes legend long 
  )

ggsave(
  filename = file.path(path, "figures/Age_Heatmap.png"),
  plot = last_plot(),
  width = 12,
  height = 22,
  dpi = 300
)

# Note: Ages 4-9: devi_anger (dev0100), devi_arg (dev0102), devi_par (dev0101) are self-reported (checked and confirmed in data documentation)
# No overlap in variables between ages 4-9 and 10+
# Small sample sizes and sparse data coverage for ages 4-10 and 27–33 


# Exclude ages 4–10 for the factor analysis, measurement discontinuity and sparse coverage 
df_ext <- df_ext %>%
  filter(age_yrs_pq >= 11, age_yrs_pq <= 26)

# Check new age range
range(df_ext$age_yrs_pq) # correct

# Check sample size
n_distinct(df_ext$pid) # updated sample size: 7460 pids


# ============================================================
# 6. OUTLIER INSPECTION AND CLEANING
# ============================================================

# Count pids per value per externalizing variable to inspect data distribution and potential outliers
df_ext_value_counts <- df_ext %>%
  select(-age_mon_pq, -age_yrs_pq, -fid) %>% # deselect demographics
  mutate(across(-pid, ~ as.numeric(.x))) %>%
  pivot_longer(cols = -pid, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  
  group_by(variable, value) %>% # count unique pids per value per variable
  summarise(n_pids = n_distinct(pid), .groups = "drop") %>%
  arrange(value) %>% # sort by value
  
  pivot_wider(names_from = value, values_from = n_pids) %>% # pivot values into columns
  arrange(variable)  # sort by variable group 

# Visual inspection
View(df_ext_value_counts)

# 1. Inspect variables' skewness
# Define value columns
val_cols <- setdiff(names(df_ext_value_counts), "variable")
val_nums <- as.numeric(val_cols)

# Identify skewed variables
skew_summary_counts <- df_ext_value_counts %>%
  rowwise() %>%
  mutate(
    w = list(replace_na(as.numeric(c_across(all_of(val_cols))), 0)), # reconstruct weights from wide frequency table format, wrap counts vector into list column
    n = sum(unlist(w)), # unwrap list back into vector
    
    w_mean = sum(val_nums * unlist(w)) / n, # weighted mean
    w_sd   = sqrt(sum(unlist(w) * (val_nums - w_mean)^2) / (n - 1)), # weighted sd
    
    skewness = (n / ((n - 1) * (n - 2))) * sum(unlist(w) * ((val_nums - w_mean) / w_sd)^3) #weighted skewness formula
  ) %>%
  select(variable, skewness) %>%
  filter(skewness > 1) %>% 
  arrange(desc(skewness))

print(n = Inf, skew_summary_counts)

# Remove outliers above Q3 + 3 * IQR for continuous skewed variables (robust to skewness because it doesn't assume normality, unlike SD-based methods; robust to zero inflation compared to log transformation + SD
# Define only continuous skewed variables to clean (exclude likert scales, yes/no etc.)
vars_skewed <- c("del_bully12", "del_shoplift12", "del_cyberb12", "del_fair_dg12", 
                 "del_steal24", "del_illeg_dl12", "del_skip_freq12", "del_vand12",
                 "del_graff12", "del_threat_rob12", "del_away_ovnt12", "drugs_freq24", 
                 "alc_n_wine_exc", "drugs_freq12", "alc_n_hp_exc", "alc_n_beer_exc")

# Identify upper fence 3 * IQR outliers for skewed variables
df_ext_value_counts_outliers_skewed <- df_ext_value_counts %>%
  filter(variable %in% vars_skewed) %>%
  rowwise() %>%
  mutate(
    w_vec = list(replace_na(as.numeric(c_across(all_of(val_cols))), 0)),
    w_q1  = Hmisc::wtd.quantile(val_nums, weights = unlist(w_vec), probs = 0.25),
    w_q3  = Hmisc::wtd.quantile(val_nums, weights = unlist(w_vec), probs = 0.75),
    w_iqr = w_q3 - w_q1,
    upper_fence = w_q3 + (3 * w_iqr),
    
    # Flags the outliers by keeping only values > fence
    across(all_of(val_cols),
           ~ if_else(!is.na(.) & as.numeric(cur_column()) > upper_fence, ., NA))
  ) %>%
  ungroup()

# Replace identified IQR outliers with NA 
df_ext <- df_ext %>%
  mutate(across(
    all_of(vars_skewed), 
    ~ {
      var_name <- cur_column()
      
      # Pull outlier value names from IQR summary
      outlier_vals <- df_ext_value_counts_outliers_skewed %>% 
        filter(variable == var_name) %>%
        select(-variable) %>% 
        select(where(~ any(!is.na(.)))) %>% 
        names() %>%
        as.numeric()
      
      # Replace outlier value with NA
      if_else(. %in% outlier_vals, NA, .) 
    }
  ))

# Check outlier removal by checking new maximum values
df_ext %>%
  select(all_of(vars_skewed)) %>%
  summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
  pivot_longer(everything()) %>%
  print(n = 20)  # worked as expected


# 2. Inspect normally distributed variables
# Get all column names from your data
all_vars <- df_ext_value_counts$variable

# Define normally distributed variables substracting the skewed variables
vars_normal <- setdiff(all_vars, vars_skewed)
vars_normal

# Identify > 4 * SD outliers for normally distributed variables
df_ext_value_counts_outliers_normal <- df_ext_value_counts %>%
  filter(variable %in% vars_normal) %>%
  rowwise() %>%
  mutate(
    w_list = list(replace_na(as.numeric(c_across(all_of(val_cols))), 0)), # wrap weights in list
    
    w_mean = sum(val_nums * unlist(w_list)) / sum(unlist(w_list)), # unlist to calculate mean
    w_sd   = sqrt(sum(unlist(w_list) * (val_nums - w_mean)^2) / (sum(unlist(w_list)) - 1)), # calculate sd
    
    across(all_of(val_cols),
           ~ if_else(!is.na(.) & abs(as.numeric(cur_column()) - w_mean) > 4 * w_sd, ., NA))  # Identify outliers > 4*SD from mean on both ends
  ) %>%
  ungroup() %>%
  select(-w_list, -w_mean, -w_sd) %>%  # drop temporary list and stats 
  select(variable, where(~ any(!is.na(.)))) %>%  # keep only rows/columns with outliers
  filter(if_any(-variable, ~ !is.na(.)))

# Visually inspect
head(df_ext_value_counts_outliers_normal, 10)

# Define continuous normal vars to clean (exclude likert scales, yes/no etc.)
vars_to_clean_normal <- c("alc_init", "alc_init_exc", "drugs_init", "smok_init"  )

# Remove normally distributed outliers > 4SD
df_ext <- df_ext %>%
  mutate(across(
    all_of(vars_to_clean_normal), 
    ~ {
      var_name <- cur_column()
      
      outlier_vals <- df_ext_value_counts_outliers_normal %>%  # find outlier values in outlier df
        filter(variable == var_name) %>%
        select(-variable) %>% 
        select(where(~ any(!is.na(.)))) %>%  # select any non-na values i.e. outlier values
        names() %>%
        as.numeric()
      
      if_else(as.numeric(.) %in% outlier_vals, NA, .)  # replace outlier values with NA
    }
  ))

# Check max/min for normal vars to confirm correct removal
df_ext %>%
  select(all_of(vars_to_clean_normal)) %>%
  summarise(across(everything(), list(
    min = ~ min(.x, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  print(n = Inf)

# Count pids with updated available data per variable
df_ext_counts_final <- df_ext %>%
  mutate(across(-pid, ~ as.numeric(.x))) %>%  # convert labelled → numeric
  pivot_longer(cols = -pid, names_to = "variable", values_to = "value") %>%  # long format
  filter(!is.na(value)) %>%  # only non missing values
  group_by(variable) %>% # count per variable
  summarise(n_unique_pids = n_distinct(pid), .groups = "drop") %>%  # count unique pids
  arrange(desc(n_unique_pids)) # sort by pids count

View(df_ext_counts_final) # completeness for age, sex, fid confirmed (7460 pids)

# Final visual check cleaned data
df_ext_value_counts_final <- df_ext %>%
  select(-age_mon_pq, -age_yrs_pq, -fid) %>% # deselect demographics
  mutate(across(-pid, ~ as.numeric(.x))) %>%
  pivot_longer(cols = -pid, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  
  group_by(variable, value) %>% # count unique pids per value per variable
  summarise(n_pids = n_distinct(pid), .groups = "drop") %>%
  arrange(value) %>% # sort by value
  
  pivot_wider(names_from = value, values_from = n_pids) %>% # pivot values into columns
  arrange(variable)  # sort by variable group 

View(df_ext_value_counts_final) # outliers successfully removed

# ============================================================
# 7. REVERSE CODING
# ============================================================

# Create a new dataframe so that high values indicate high externalizing for variables with good coverage across age
df_ext <- df_ext %>%
  mutate(
    # attention (flip scale 1-3 : 1=Low -> 3=High)
    adhd_act   = 4 - adhd_act,
    adhd_tasks = 4 - adhd_tasks,
    
    # delinquency, drugs (flip 1=Yes/2=No to 1=No/2=Yes)
    del_illeg_dl_ever   = 3 - del_illeg_dl_ever,
    del_shoplift_ever   = 3 - del_shoplift_ever,
    del_cyberb_ever   = 3 - del_cyberb_ever,
    del_vand_ever    = 3 - del_vand_ever,
    del_drink_drive_ever = 3 - del_drink_drive_ever,
    drugs_ever = 3 - drugs_ever,   
    
    del_compl = 4 - del_compl, # 3=applies completely/1=doesn't apply at all to # 1=applies completely/3=doesn't apply at all
    
    # smok, alc (flip 1=Freq -> 6=Never)
    smok = 7 - smok,      
    alc_freq_exc  = 7 - alc_freq_exc,   
    
    # alc_init, smok_init
    alc_init = 26 - alc_init,    # Lower age of first drink = Higher score (33 max age)
    smok_init = 26 - smok_init
  )

# ============================================================
# 8. DESCRIPTIVE STATISTICS
# ============================================================

# Get descriptives for all participants across all timepoints
df_descriptives <- df_ext %>%
  select(pid,age_mon_pq, age_yrs_pq, sex, zygosity) %>% # select variables for descriptives
  mutate(across(everything(), ~ if_else(. %in% missing_codes, NA, as.numeric(.)))) %>%
  summarise(
    across(c(age_mon_pq, age_yrs_pq), list( # get descriptives for age in months and years
      n      = ~ sum(!is.na(.)),
      mean   = ~ mean(., na.rm = TRUE),
      sd     = ~ sd(., na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      min    = ~ min(., na.rm = TRUE),
      max    = ~ max(., na.rm = TRUE)
    )),
    
    # pid-level variables: compute per pid not per row (multiple timpoints per pid)
    sex_n          = n_distinct(pid[!is.na(sex)]),
    sex_pct_female = mean(sex[!duplicated(pid)] == 2, na.rm = TRUE) * 100, # percentage of female pids
    zyg_n          = n_distinct(pid[!is.na(zygosity)]),
    zyg_pct_mz     = mean(zygosity[!duplicated(pid)] == 1, na.rm = TRUE) * 100 # percentage of monozygotic pids
    
  )

View(df_descriptives)



# ============================================================
# 9. SAVE DATA
# ============================================================

# Export descriptive statistics
write.csv(df_descriptives, file.path(path, "02_twinLife_externalizing_descriptives.csv"), row.names = FALSE)


# Export cleaned externalizing dataset
write.csv(df_ext, file.path(path, "02_twinLife_externalizing_data_cleaned.csv"), row.names = FALSE)
save(df_ext, file = file.path(path, "02_twinLife_externalizing_data_cleaned.rda"))




