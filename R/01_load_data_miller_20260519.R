########################################################################
#### Externalizing — Phenotypic data mining in TwinLife             #### 
#### 25.02.2026 - 19.05.2026                                        ####
#### Load TwinLife demographic, zygosity and externalizing data     ####
#### main Author: Anna Miller                                       #### 
########################################################################

# ============================================================
# 0. IMPORT LIBRARIES
# ============================================================

library(haven) # to read dta file provided by TwinLife
library(dplyr) # filter and select data

# ============================================================
# 1. DATA LOADING
# ============================================================

# Define path
path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2")

# Load Twinlife dataset provided by TwinLife
df_twinlife <- read_dta(file.path(path, "20260320_TwinLife_Externalizing_Phenotype.dta"))


# ============================================================
# 2. FILTER SELF-REPORTED TWIN DATA
# ============================================================

# Define respondents variables according to data documentation
respondents <- c(
  "pid",
  "ptyp_1",
  "ptyp_2",
  "ptyp_3",
  "ptyp_4",
  "ptyp_5",
  "ptyp_6",
  "ptyp_7",
  "ptyp_8",
  "ptyp_9"
)

# Create df to inspect
df_respondents <- df_twinlife[, respondents]
head(df_respondents)

df_respondents_sorted <- df_respondents[order(df_respondents$pid), ]
head(df_respondents_sorted, 20) # multiple rows per pid, not all contain ptyp data

# Count number of pids for first respondent type variable to confirm 
# data availability for 4096 twin pairs as expected from data documentation
df_ptyp_counts <- df_twinlife %>%
  filter(!is.na(ptyp_1)) %>%                                   
  group_by(ptyp_1) %>%                                        
  summarise(pid_count = n_distinct(pid))    

df_ptyp_counts # data availability for 4096 twin pairs (ptyp 1 and 2) confirmed

# Check for consistency across respondent type variables
df_mismatch <- df_respondents %>%
  mutate(
    mismatch = {
      m <- as.matrix(select(., starts_with("ptyp_")))  
      m[m == -95] <- NA                                # set missing code -95 as NA
      m[m == 201] <- 200                               # merge respondent type 201 (sibling not participating) with 200 (sibling)
      rowSums(m != m[, 1], na.rm = TRUE) > 0           # find mismatches
    }
  )

df_mismatch %>%
  filter(mismatch) # No mismatch in respondent type twin (1, 2)

# Check whether all pids have at least one non-NA ptyp_1 value, 
# or if other ptyp variables must be used to fill gaps
pids_all_ptyp_1_na <- df_respondents %>%
  group_by(pid) %>%
  summarise(all_na = all(is.na(ptyp_1))) %>%  # find pids for which no row contains ptyp_1
  filter(all_na) %>%                          
  pull(pid) %>%
  as.numeric()

pids_all_ptyp_1_na # all pids have at least one ptyp_1 value

# Create df with pids and ptyp_1 to check for duplicates
df_ptyp_1 <- df_respondents %>%
  filter(!is.na(ptyp_1) & !(ptyp_1 != -95)) %>%  
  select(pid, ptyp_1)                                        

head(df_ptyp_1)

# Check for duplicate pids to verify consistency
df_ptyp_1 %>%
  group_by(pid) %>%
  filter(n() > 1)   # no duplicates

# Filter df_twinlife for twin respondents to filter for self-reported data
df_twin_pids <- df_twinlife %>%
  filter(ptyp_1 %in% c(1, 2))

# Create df for twin respondents with self-reported data only
df_twins <- df_twinlife %>%
  filter(pid %in% df_twin_pids$pid)

# Verify complete number of pids per twin respondent type
df_ptyp_counts <- df_twins %>%
  filter(!is.na(ptyp_1)) %>%                                
  group_by(ptyp_1) %>%                                      
  summarise(pid_count = n_distinct(pid), .groups = "drop")  

df_ptyp_counts # dataset complete with 4096 twins


# ============================================================
# 3. ZYGOSITY DATA
# ============================================================

# Define zygosity variables to find best variable to use
zyg <- c(
  "pid", #person ID
  "zyg0100", #zygosity: result questionnaire in F2F1 (gen)
  "zyg0101", #zygosity: result DNA in F2F1 (gen)
  "zyg0102", #zygosity: result questionnaire and DNA in F2F1 (gen)
  "zyg0109", #value of zygosity discriminant function (gen)
  "zyg0111", #zygosity: result saliva sample in F2F3 (gen)
  "zyg0112", #zygosity: result questionnaire and DNA in F2F1 & F2F3 (gen)
  "zyg0700", #assessment: zygosity of twins (m/f)
  "zyg0700_t", #assessment twin 1: zygosity of twins (gen)
  "zyg0700_u", #assessment twin 2: zygosity of twins (gen)
  "sdc0100" #subsequently completed zygosity questionnaire (gen)
 )

# Create df with zygosity variables for all twin respondents to check zyg data availablity
df_zyg <- df_twins[, zyg]
head(df_zyg) # questionnaire and DNA/saliva results have different values for some pids

# Check which values variables have to check completeness of missing codes defined above
unique(df_zyg$zyg0101)
unique(df_zyg$zyg0111)

# Define missing codes
missing_codes <- c(-99, -98, -95, -94, -93, -92, -90, -89, -87)

# Replace missing codes with NA
df_zyg <- df_zyg %>%
  mutate(across(
    everything(),
    ~ replace(.x, as.numeric(.x) %in% missing_codes, NA)
  ))

# Check whether zyg0101 (DNA F2F1) and zyg0111 (saliva F2F3) mismatch when not missing
df_mismatch_zyg <- df_zyg %>%
  mutate(
    mismatch_zyg = {
      dna   <- ifelse(is.na(zyg0101), NA, zyg0101)
      sal   <- ifelse(is.na(zyg0111), NA, zyg0111)
      !is.na(dna) & !is.na(sal) & dna != sal
    }
  )

# Filter for mismatches to inspect
df_mismatch_zyg <- df_mismatch_zyg %>% 
  filter(mismatch_zyg) %>% 
  select(pid, zyg0101, zyg0111, mismatch_zyg) 

df_mismatch_zyg # mismatch for 2 pids

# Check all zyg variables for mismatching pids
df_check <- df_zyg %>% filter(pid %in% df_mismatch_zyg$pid)
head(df_check, 20)

# Conservatively remove zygosity data due to inconsistency for these 2 pids
df_zyg <- df_zyg %>%
  mutate(
    across(
      -pid,
      ~ replace(.x, pid %in% df_mismatch_zyg$pid, NA)
    )
  )

# Check if zyg0111 is not available but zyg0101 is for any pid
df_zyg %>%
  group_by(pid) %>%
  summarise(
    has_zyg0101 = any(!is.na(zyg0101)),
    has_zyg0111 = any(!is.na(zyg0111))
  ) %>%
  filter(has_zyg0101 & !has_zyg0111) %>%
  nrow() # 272 without zyg0111 but with zyg0101
  
# Since also zyg0101 is not available for all pids: 
# check if latest questionnaire (zyg0112) is available for all twins
df_check <- df_zyg %>%
  group_by(pid) %>%
  filter(!any(!is.na(zyg0112))) %>%
  ungroup()                            

head(df_check, 20) 
n_distinct(df_check$pid) # None of all the zygosity variables available for 8 pids 

# Create df with final zyg variable for all twin pids
df_zyg_final <- df_zyg %>%
  group_by(pid) %>%                               # zyg variables per pid
  summarise(
    across(c(zyg0111, zyg0101, zyg0112), ~ {      # check availability zyg vars
      x <- .[!is.na(.)]
      if (length(x) > 0) last(x) else NA
    })
  ) %>%
  mutate(zygosity = coalesce(zyg0111, zyg0101, zyg0112)) %>% # use latest zyg0111 (saliva) if available, 
                                                             # otherwise zyg0101 (dna) if available,
                                                             # otherwise zyg0112 (latest questionnare)
  filter(pid %in% df_twins$pid) %>%
  select(pid, zygosity)

# Check if there are any pids without zygosity data
df_zyg_final %>%
  filter(is.na(zygosity)) %>%
  nrow() # 8 pids without (consistent) zygosity data will have to be excluded from twin model


# ============================================================
# 4. EXTERNALIZING DATA & DEMOGRAPHICS
# ============================================================

# Define variables for demographics and externalizing related self-reported measures according to data documentation
ext <- c(
  "pid", # person ID
  "fid", # family ID
  "wid", # data collection ID - consecutively numbered (gen)
  "sex", #  self-reported sex
  
  # Birth cohort, age and years of data collection
  "cgr", # twin birth cohort
  "fpr0104", # month of birth
  "fpr0105", # year of birth
  
  "age0100", # age in years on the date of the family questionnaire (gen)
  "age0101", # age in months "
  "age0100t", # age twins in years on the date of the family questionnaire (gen),
  "yea_pq", # year of individual questionnaire
  "mon_pq", # month of individual questionnaire
  "yea_fq", # year of family questionnaire
  "mon_fq", # month of family questionnaire
  
  # Smoking
  "del1200", # Ever smoked
  "hbe0100", # Smoking behavior 
  "del1201", # Age of smoking initiation
  "hbe0101", # Age of smoking initiation
  "del1203", # Smoking frequency
  "del1204", # Smoking alone/ in group
  "hbe0102", # Smoking alone/ in group
  
  # Alcohol
  "del1300", # Ever consumed alcohol
  "del1304", # Ever excessively 
  "del1301", # Age of drinking initiation
  "hbe0200", # Age of drinking initiation
  "del1305", # Age of first time being drunk
  "del1309", # Frequency excessive alcohol consumption
  "hbe0202", # Frequency excessive alcohol consumption
  "del1306", # Frequency excessive alcohol consumption: Number of beverages - beer 0.3
  "hbe0210", # Frequency excessive alcohol consumption: Number of glasses - beer
  "hbe0220", # Frequency excessive alcohol consumption: Number of glasses - wine
  "hbe0230", # Frequency excessive alcohol consumption: Number of glasses - high proof alcohol
  "dia1000", # Alcohol addiction diagnosed
  
  # Drug consumption
  "del1500", # Ever consumed drugs
  "del1501", # Frequency drug consumption last 12 months
  "del1505", # Frequency drug consumption last 24 months
  "del1502", # Age of first drug consumption
  "del1503", # Consuming drugs alone/ in group
  
  # Risk taking
  "per0200", # Willingness to take risks 
  "ris0100", # Risk tolerance
  
  # Self-regulation (BISS Scale)
  "srg0100", # New ideas distract
  "srg0200", # Pursue different goals
  "srg0300", # Changing interests
  "srg0400", # Do bad fun things
  "srg0500", # Fun distracts from job
  "srg0600", # Wish for more self discipline
  
  # Hyperactivity / Inattention (SDQ)
  "ext0100", # Restlessness
  "ext0101", # Moving, fidgeting
  "ext0102", # Distracted, unfocused
  "ext0103", # Thinks before acting
  "ext0104", # Finishes tasks, able to concentrate
  
  # ADHD/ ADD diagnosed last 12 months
  "dia4700", 
  
  # Deviant behavior
  "dev0100", # Fits of anger
  "ext0105", # Fits of anger
  "dev0102", # Listening to parents
  "dev0103", # Arguing with others, bullying
  
  # Aggression
  "ext0107", # attack other physically
  
  # Rule-Breaking/ Delinquency
  "del1700", # Drive drunk or under drugs ever 
  "del1704", # Drive drunk or under drugs last 24 months
  "ext0109", # Stealing
  "del0500", # Stealing ever
  "del0501", # Stealing last 24 months
  "del0502", # Stealing frequency
  "del0100", # fair dodging
  "del0101", # fair dodging last 12 months frequency
  "del0200", # away over night ever
  "del0201", # away over night last 12 months frequency
  "del0202", # away over night alone or in group
  "del0300", # skip school ever
  "del0301", # skip school last 12 months frequency
  "del0302", # skip school alone or in group
  "del0303", # skip school several days
  "del0400", # shoplift ever
  "del0401", # shoplift last 12 months frequency
  "del0402", # shoplift alone or in group
  "del0600", # graffiti ever
  "del0601", # graffiti last 12 months frequency
  "del0602", # graffiti alone or in group
  "del0700", # vandalism ever
  "del0701", # vandalism last 12 months frequency
  "del0702", # vandalism alone or in group
  "del1000", # bullying ever
  "del1001", # bullying last 12 months frequency
  "del1002", # bullying alone or in group
  "del0800", # cyber bullying ever
  "del0801", # cyber bullying last 12 months frequency
  "del0900", # illegal downloading ever
  "del0901", # illegal downloading last 12 months frequency
  "del1100", # threatening/ robbing ever
  "del1101", # threatening/ robbing last 12 months frequency
  "del1102", # threatening/ robbing alone or in group
  "ext0106", # Compliance, obedience
  "ext0108", # Lying
  
  # Sexual behavior
  "seo0900" # STD protection
)  

# Create externalizing df with self-reported twins data
df_ext <- df_twins[, ext]
head(df_ext)

# Join with zygosity data
df_ext <- df_ext %>%
  left_join(df_zyg_final %>% select(pid, zygosity), by = "pid")

# Check completeness of join
df_ext %>%
  group_by(pid) %>%
  filter(!any(!is.na(zygosity))) %>%     # filter for pids without zygosity data
  distinct(pid) %>%
  nrow()                                  # complete join, 8 pids without zygosity data

# Remove pids without zygosity data
df_ext <- df_ext %>%
  group_by(pid) %>%
  filter(any(!is.na(zygosity))) %>%
  ungroup()
  
# Rename variables
df_ext <- df_ext %>%
  rename(
    # Birth and age
    birth_month = fpr0104,
    birth_year = fpr0105,
    age_yrs_fq = age0100, 
    age_mon_fq = age0101,
    age_twins_yrs_fq = age0100t,
    
    # Smoking
    smok = hbe0100,
    smok_ever = del1200, 
    smok_init = del1201,
    smok_init2 = hbe0101,
    smok_freq = del1203,
    smok_solo = del1204,
    smok_solo2 = hbe0102,
    
    # Alcohol
    alc_ever = del1300,
    alc_ever_exc = del1304,
    alc_init = del1301,
    alc_init2 = hbe0200,
    alc_init_exc = del1305,
    alc_freq_exc = del1309,
    alc_freq_exc2 = hbe0202,
    alc_n_beer_exc = del1306, # Number of beers 0.3
    alc_n_beer_exc2 = hbe0210, # Number of glasses - beer
    alc_n_wine_exc = hbe0220,
    alc_n_hp_exc = hbe0230, 
    alc_diag = dia1000, 
    
    # Drug consumption
    
    drugs_ever = del1500,
    drugs_freq12 = del1501, 
    drugs_freq24 = del1505, 
    drugs_init = del1502,
    drugs_solo = del1503, 
    
    # Risk taking
    risk_will = per0200,
    risk_tol = ris0100,
    
    # Hyperactivity / Inattention (SDQ) 
    adhd_restless = ext0100,
    adhd_move = ext0101,
    adhd_unfoc = ext0102,
    adhd_act = ext0103,
    adhd_tasks = ext0104,
    
    # ADHD/ ADD diagnosed last 12 months
    adhd_add_diag = dia4700,
    
    # Deviant behavior
    devi_anger = dev0100, 
    devi_anger2 = ext0105, 
    devi_par = dev0102, 
    devi_arg = dev0103, 
    
    # Aggression
    agg = ext0107, # attack other physically
    
    # Rule-Breaking/ Delinquency
    del_drink_drive_ever = del1700, 
    del_drink_drive24 = del1704, 
    del_steal = ext0109, 
    del_steal_ever = del0500,
    del_steal24 = del0501, 
    del_steal_freq = del0502, 
    del_fair_dg_ever = del0100, 
    del_fair_dg12 = del0101, 
    del_away_ovnt = del0200, 
    del_away_ovnt12 = del0201, 
    del_away_ovnt_solo = del0202,
    del_skip_ever = del0300, 
    del_skip_freq12 = del0301, 
    del_skip_solo = del0302, 
    del_skip_days = del0303,
    del_shoplift_ever = del0400, 
    del_shoplift12 = del0401, 
    del_shoplift_solo = del0402, 
    del_graff_ever = del0600, 
    del_graff12 = del0601, 
    del_graff_solo = del0602, 
    del_vand_ever = del0700, 
    del_vand12 = del0701, 
    del_vand_solo = del0702, 
    del_bully_ever = del1000,
    del_bully12 = del1001,
    del_bully_solo = del1002,
    del_cyberb_ever = del0800, 
    del_cyberb12 = del0801, 
    del_illeg_dl_ever = del0900, 
    del_illeg_dl12 = del0901, 
    del_threat_rob_ever = del1100, 
    del_threat_rob12 = del1101,
    del_threat_rob_solo = del1102, 
    del_compl = ext0106,
    del_lie = ext0108,
    
    # Sexual behavior
    sex_prot = seo0900
  )


# ============================================================
# 5. SAVE DATA
# ============================================================

# Export twins self-reported externalizing dataset
write.csv(df_ext, file.path(path, "Processed_Data", "01_twinLife_externalizing_data.csv"), row.names = FALSE)
save(df_ext, file = file.path(path, "Processed_Data/01_twinLife_externalizing_data.rda"))
