########################################################################
#### Externalizing — Timepoint Selection in TwinLife                #### 
#### 22.03.2026 - 19.05.2026                                        ####
#### Cohort-based cross-sectional sampling                          ####
#### main Author: Anna Miller                                       #### 
########################################################################

# ============================================================
# 0. IMPORT LIBRARY FOR DATA MANIPULATION
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 1. LOAD DATA
# ============================================================

path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data/")

# Load preprocessed data
load(file.path(path, "02_twinLife_externalizing_data_cleaned.rda"))

n_distinct(df_ext$pid) # confirm correct sample size 7460 pids
colnames(df_ext)


# ============================================================
# 2. PREPARE DATA
# ============================================================

# Define demographic variables 
demo_vars <- c("pid", "fid", "sex", "zygosity", "age_mon_pq", "age_yrs_pq", "yea_pq", "mon_pq")

# Remove rows that have none of the externalizing related variables but only demographic
df_ext_filtered <- df_ext %>%
  filter(if_any(-all_of(demo_vars), ~ !is.na(.)))

n_distinct(df_ext_filtered$pid) # sample size 7403 pids


# Srg available in some rows where none of the externalizing related variables are
# Define srg vars
srg_vars <- c("srg0100", "srg0200", "srg0300", 
              "srg0400", "srg0500", "srg0600")

# Define others
other_vars <- setdiff(names(df_ext), c("pid", "age_yrs_pq", srg_vars))

# Check if there are rows for same pids at same age with conflicting srg
df_ext_filtered %>%
  group_by(pid, age_yrs_pq) %>%
  summarise(
    across(everything(), ~ n_distinct(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -c(pid, age_yrs_pq),
    names_to = "variable",
    values_to = "n_values"
  ) %>%
  filter(n_values > 1) %>%
  count(variable, sort = TRUE)
  # Conflicting data only for srg and age_mon_pq

# Merge rows per pid per age and take mean of different srg values at same age
df_ext_merged <- df_ext_filtered %>%
  group_by(pid, age_yrs_pq) %>%
  summarise(
    across(all_of(other_vars), ~ {
      x <- na.omit(.x)
      if (length(x) == 0) NA else x[1]
    }),
    across(all_of(srg_vars), ~ {
      x <- na.omit(as.numeric(.x))
      if (length(x) == 0) NA_real_ else mean(x)
    }),
    .groups = "drop"
  )

# Check merged
df_ext_merged %>%
  count(pid, age_yrs_pq) %>%
  filter(n > 1) # no more conflicts


# Count available observations per age year across all pids
age_counts <- df_ext_merged %>%
  count(age_yrs_pq, name = "n_pids")

age_counts


# ============================================================
# 3. VISUALIZATIONS OF DATA COLLECTION ACROSS TIME
# ============================================================

# 1. Visualize repeated data collection
# Prepare timepoints to plot data collection across years (months not available in all rows)
plot_repeated_years <- df_ext_merged %>%
  distinct(pid, yea_pq) %>%
  count(pid, name = "n_years") %>%
  count(n_years, name = "n_participants")

# Plot number of participants with their number of years they participated in
plot_repeated_data <- ggplot(
  plot_repeated_years,
  aes(x = n_years, y = n_participants)
) +
  geom_col(fill = "#7FBFA2") +
  scale_x_continuous(
    breaks = sort(unique(plot_repeated_years$n_years))
  ) +
  labs(
    title = "Repeated TwinLife data-collection across years",
    x = "Number of years with available data per participant",
    y = "Number of participants"
  ) +
  theme_minimal(base_size = 16) # minimal theme with increased text size

plot_repeated_data

# sanity check sum of number of participants
sum(plot_repeated_years$n_participants) # correct

# Save
ggsave(
  filename = file.path(path, "figures/Repeated_data_collection.png"),
  plot = plot_repeated_data,
  width = 12,
  height = 8,
  dpi = 300
)

# 2. Visualize data collection by years
# Count unique participants data was collected for per year
plot_year_counts <- df_ext_merged %>%
  distinct(pid, yea_pq) %>%
  count(yea_pq, name = "n_participants")

# Plot data collection by year
plot_collection_years <- ggplot(
  plot_year_counts,
  aes(x = factor(yea_pq), y = n_participants)
) +
  geom_col(fill = "#3845AD") +
  geom_text(
    aes(label = n_participants),
    vjust = -0.4,
    size = 5
  ) +
  labs(
    title = "TwinLife data collection across years",
    x = "Year of data collection",
    y = "Number of participants"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  expand_limits(
    y = max(plot_year_counts$n_participants) * 1.08
  )

plot_collection_years

# Save
ggsave(
  filename = file.path(path, "figures/Data_collection_years.png"),
  plot = plot_collection_years,
  width = 12,
  height = 8,
  dpi = 300
)

# ============================================================
# 4. CROSS-SECTIONAL SAMPLING
# ============================================================

# Count available observations per age year across all twin pairs
twins_age_counts <- df_ext_merged %>%
  group_by(fid, pid) %>%
  summarise(ages = list(unique(age_yrs_pq)), .groups = "drop") %>%
  group_by(fid) %>%
  summarise(age = list(Reduce(intersect, ages)), .groups = "drop") %>%
  tidyr::unnest(age) %>%
  count(age, name = "n_twin_pairs")

twins_age_counts

# Check if there are twin pairs with no common age
twins_no_common_age <- df_ext_merged %>%
  group_by(fid, pid) %>%
  summarise(ages = list(unique(age_yrs_pq)), .groups = "drop") %>%
  group_by(fid) %>%
  summarise(
    common_ages = list(Reduce(intersect, ages)),
    .groups = "drop"
  ) %>%
  filter(lengths(common_ages) == 0) %>%
  select(fid)

twins_no_common_age
# 11 twins, will need else clause to select different timepoints for twins with no common timepoint
  
# Set median count as target for more even distribution
target <- df_ext_merged %>%
  count(age_yrs_pq) %>%
  pull(n) %>%
  median()

# Initialize tracker that counts age
age_tracker <- df_ext_merged %>%
  distinct(age_yrs_pq) %>%
  mutate(n_selected = 0)

# Initialize selected rows
selected_rows <- tibble()

# Randomly sample fids to sample per twin pair
set.seed(123) # set seed for reproducability
fids <- sample(unique(df_ext_merged$fid))

# Loop through randomized twin pairs (fids) to select timepoints at ages below target
for (f in fids) {
  fid_data <- df_ext_merged %>%
    filter(fid == f) %>%
    left_join(age_tracker, by = "age_yrs_pq")
  
  # find ages available for all pids in this family
  common_ages <- Reduce(intersect, split(fid_data$age_yrs_pq, fid_data$pid))
  
  if (length(common_ages) > 0) {
    # select best common age based on tracker
    candidates <- fid_data %>% filter(age_yrs_pq %in% common_ages)
    not_at_target <- candidates %>% filter(n_selected < target)
    
    best_age <- if (nrow(not_at_target) > 0) {
      not_at_target %>% slice_min(order_by = n_selected, n = 1, with_ties = FALSE) %>% pull(age_yrs_pq)
    } else {
      candidates %>% slice_min(order_by = n_selected, n = 1, with_ties = FALSE) %>% pull(age_yrs_pq)
    }
    
    selected <- fid_data %>% filter(age_yrs_pq == best_age) %>% select(-n_selected)
    
  } else {
    # if no common age available then selection per pid i.e. per twin individual
    selected <- fid_data %>%
      group_by(pid) %>%
      group_modify(~ {
        not_at_target <- .x %>% filter(n_selected < target)
        if (nrow(not_at_target) > 0) {
          not_at_target %>% slice_min(order_by = n_selected, n = 1, with_ties = FALSE)
        } else {
          .x %>% slice_min(order_by = n_selected, n = 1, with_ties = FALSE)
        }
      }) %>%
      ungroup() %>%
      select(-n_selected)
  }
  
  # update age tracker based on selected data
  age_tracker <- age_tracker %>%
    mutate(n_selected = n_selected + sapply(age_yrs_pq, function(a) sum(selected$age_yrs_pq == a)))
  
  selected_rows <- bind_rows(selected_rows, selected)
}


# Create new dataset with selected timepoints
df_ext_cross <- selected_rows

# Check age distribution
twins_age_counts_new <- df_ext_cross %>%
  count(age_yrs_pq) %>%
  print(n = Inf)

# Compare age distributions before after
twins_age_counts_merged <- twins_age_counts_new %>%
  left_join(twins_age_counts, by = c("age_yrs_pq" = "age")) %>%
  rename(                       
    n_twin_pairs_after = n,
    n_twin_pairs_before = n_twin_pairs
  ) %>%
  relocate(n_twin_pairs_before, n_twin_pairs_after, .after = age_yrs_pq) # change order

twins_age_counts_merged

# Check completeness in timepoint selection
n_distinct(df_ext_filtered$pid)
n_distinct(df_ext_cross$pid) # all pids represented
nrow(df_ext_cross) # one row (timepoint) per pid 

# Check if same age for twin pairs except for 11 twins with no common age
df_check_no_common_age <- df_ext_cross %>%
  group_by(fid) %>%
  filter(n() == 2) %>%
  summarise(
    n_age_yrs = n_distinct(age_yrs_pq),
    .groups = "drop"
  ) %>%
  filter(n_age_yrs > 1 )

df_check_no_common_age # correctly selected same age except for 11 twins with no common age

# ============================================================
# 5. SAVE DATA
# ============================================================

# Export externalizing cross-sectional sampling dataset
write.csv(df_ext_cross, file.path(path, "03_twinlife_externalizing_cross_sampled.csv"), row.names = FALSE)
save(df_ext_cross, file = file.path(path, "03_twinlife_externalizing_cross_sampled.rda"))
