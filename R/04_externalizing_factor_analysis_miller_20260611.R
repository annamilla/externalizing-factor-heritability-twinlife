########################################################################
#### Externalizing — Factor Analysis and Invariance in TwinLife     #### 
#### 01.04.2026 - 08.05.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

library(dplyr)
library(lavaan)
library(tidyr)
library(ggplot2)
library(corrplot)
library(semTools)
library(semPlot)


# ============================================================
# 1. DATA LOADING 
# ============================================================

# Define path for data files
path <- ("/Volumes/1000-twinlife/private/data/2026_TwinLife_Externalizing/20260320_TwinLife_Externalizing_Data_Transfer_V2/Processed_Data")

# Define path for factor models
path_model <- ("/Volumes/MPRG-Biosocial/Projects/04_data_analysis/012_SHIP_BASE_TwinLife_SOEP/Externalizing_Genomics/models")

# Load preprocessed data
load(file.path(path, "03_twinlife_externalizing_cross_sampled.rda"))

head(df_ext_cross)

# check completeness of columns with ext vars
colnames(df_ext_cross) 


# ============================================================
# 2. CONFIRMATORY FACTOR ANALYSIS (CFA) PREP
# ============================================================

# Define demo vars not included as indicators for CFA
demo_vars <- c("pid", "fid", "sex", "zygosity", "age_yrs_pq", "age_mon_pq")

# Define all indicators with good coverage (see first script) across age 11-26
indicators <- c(
  "adhd_act", "adhd_move", "adhd_restless", "adhd_tasks", "adhd_unfoc",
  "agg",
  "alc_freq_exc", "alc_init", "alc_init_exc", "alc_n_beer_exc",
  "alc_n_hp_exc", "alc_n_wine_exc",
  "del_compl", "del_lie", "del_steal", "del_cyberb_ever", "del_drink_drive24",
  "devi_anger",
  "drugs_ever",
  "risk_will",
  "smok", "smok_init",
  "srg0100", "srg0200", "srg0300", "srg0400", "srg0500", "srg0600")

# Check how many response categories each variable has to choose appropriate estimator (MLR/WLSMV)
df_ext_cross %>%
  select(indicators) %>%                                     # check indicators only without demographic vars
  summarise(across(everything(), ~ length(unique(na.omit(.))))) %>%  # count unique values for each indicator
  tidyr::pivot_longer(everything(),                                  # convert format
                      names_to = "variable", 
                      values_to = "n_categories") %>%
  arrange(n_categories) %>%
  print(n = Inf) # Checked variables with only 2,3 number of categories in data documentation: correct number of categories
                  # Include variables with 3+ categories in CFA only! Estimators can't handle both binary and ordinal well

# convert all variables to numeric
df_model <- df_ext_cross %>%
  # Force all columns to be numeric (treat categorical as ordinal)
  mutate(across(all_of(indicators), ~as.numeric(as.character(.))))

# Check final sample size
nrow(df_model) # --> final sample size: 7403

# Sanity check age range
range(df_ext_cross$age_yrs_pq, na.rm = TRUE) # correctly filtered for 11-26


# ============================================================
# 3. CFA MODEL SELECTION
# ============================================================

# Define helper function to fit and compare multiple models
fit_ext <- function(model, data = df_model) {
  cfa(
    model = model,
    data = data,
    estimator = "MLR",   # robust maximum likelihood 
    missing = "FIML",    # full information max likelihood to bridge gaps in data
    std.lv = TRUE,       # fix scale of latent variables while keeping all loadings free
    cluster = "fid"      # family nesting
  )
}

# 0. Load model which includes all indicators with good coverage and 3+ categories
source(file.path(path_model, "model0.R"))

# Fit model
fit_ext0 <- fit_ext(model0)
summary(fit_ext0, fit.measures = TRUE, standardized = TRUE)

# Identify variables with remaining weak coverage
coverage <- lavInspect(fit_ext0, "coverage")
round(coverage, 2)
which(coverage < .10, arr.ind = TRUE) # remove: alc_init_exc, srg0300, del_compl, alc_n_wine_exc

# Check modification indices for residual covariances and potential method effects
mi <- modificationIndices(fit_ext0, sort. = TRUE)

subset(mi, mi > 10 ) # strong residual cov adhd_restless ~~ adhd_move -> checked documentation: similar wording/content, suggests method effect
                    # also adhd_unfoc ~~ adhd_tasks -> checked documentation: similar wording/content 
                    # latent vars dis ~~ del suggest strong overlap -> model as one latent factor
                    # srg0100 ~~ srg0200 -> checked documentation: overlap in wording/content


# 1. Load redefined model1: removed variables with weak coverage and added 
# residual covariances to account for method effects
source(file.path(path_model, "model1.R"))

# Fit model
fit_ext1 <- fit_ext(model1)
summary(fit_ext1, fit.measures = TRUE, standardized = TRUE) # problem with two-indicator factor hyp and added residual covariance


# 2. Fit second model without residual cov restless/move
source(file.path(path_model, "model2.R"))

fit_ext2 <- fit_ext(model2)
summary(fit_ext2, fit.measures = TRUE, standardized = TRUE)

# Check modification indices 
mi <- modificationIndices(fit_ext2, sort. = TRUE)
subset(mi, mi > 10 ) # moderate residual cov srg0400 with alc_freq_exc and smok
                     # srg0400: doing bad things, even if they are fun -> content related to substance use!


# 3. Sensitivity model including srg with alc_freq_excand smok residual covariance
source(file.path(path_model, "model3.R"))

fit_ext3 <- fit_ext(model3)
summary(fit_ext3, fit.measures = TRUE, standardized = TRUE) # improves model fit, but robust indices fail

# compare higher order loadings to check if they are stable with and without residual cov
parameterEstimates(fit_ext3, standardized = TRUE) |>
  subset(op == "=~" & lhs == "EXT")

# compare sub and srg loadings
parameterEstimates(fit_ext3, standardized = TRUE) |>
  subset(op == "=~" & lhs %in% c("srg", "sub"))
  # improves fit but does not change higher order loadings or srg and sub loadings substantially
  # -> higher-order factor is stable also without adding this residual covariance

# 4. Sensitivity check whether hyp with only two indicators is driving the higher-order model
# Create composite mean score instead of latent hyp factor to check robustness
df_model$hyp_comp <- rowMeans(df_model[, c("adhd_restless", "adhd_move")], na.rm = TRUE) 

source(file.path(path_model, "model4.R"))

fit_ext4 <- fit_ext(model4)
summary(fit_ext4, fit.measures = TRUE, standardized = TRUE) # replacing hyp with a composite makes fit worse
                  # hyp_comp still loads on EXT , so EXT is not artifact of two indicators
                  # also core EXT structure is stable
                  # sub is weakest part of the model


# 5. Fit final model without substance use
source(file.path(path_model, "model_final.R"))

fit_ext_final <- fit_ext(model_final)
summary(fit_ext_final, fit.measures = TRUE, standardized = TRUE)

parameterEstimates(fit_ext_final, standardized = TRUE) |>
  subset(op == "=~")
  # model without substance looks cleaner
  # higher-order EXT structure stable -> EXT factor does not depend on substance use for its definition
  # con factor is still weak on agg and del_steal, but not as problematic as sub, because con loads well on EXT


# 5.2. Sensitivity model to test whether substance use is related to EXT
source(file.path(path_model, "model5.R"))

fit_ext5 <- fit_ext(model5)
summary(fit_ext5, fit.measures = TRUE, standardized = TRUE)
  # substance use is significantly associated with EXT, but moderately
  # so it's related but not strong enough to be a core component of EXT

# 5.3. Sensitivity model non-hierarchical: one EXT factor for all indicators
source(file.path(path_model, "model6.R"))

fit_ext6 <- fit_ext(model6)
summary(fit_ext6, fit.measures = TRUE, standardized = TRUE)


# Sanity check rows/pids without variables used in final model
final_model_vars <- c(
  "adhd_restless", "adhd_move", "adhd_unfoc",
  "adhd_act", "adhd_tasks",
  "srg0100", "srg0200", "srg0400", "srg0500", "srg0600",
  "agg", "devi_anger", "del_lie", "del_steal"
)

# Define rows with none of the variables
nonrelevant_vars_rows <- rowSums(!is.na(df_model[, final_model_vars])) == 0

sum(nonrelevant_vars_rows)
df_nonrelevant <- df_model[nonrelevant_vars_rows, ]
  # fine just rows with few variables such as risk_tol not used in factor analysis

# Create table fit indices from final model to report
fit_stats <- fitMeasures(fit_ext_final, c("cfi", "tli", "cfi.robust", "tli.robust", "rmsea", "rmsea.robust", "srmr"))

fit_table <- data.frame(
  Index = names(fit_stats),
  Value = round(as.numeric(fit_stats), 3)
)
fit_table


# ============================================================
# 4. FACTOR SCORES AND CORRELATIONS
# ============================================================

# Extract latent factor scores
fs <- lavPredict(fit_ext_final, 
                 type = "lv", 
                 method = "regression") # regression to calculate scores with missing data 

# Create df with variables and EXT score to compute correlations
df_corr <- df_model %>%
  mutate(EXT_score = as.numeric(fs[, "EXT"]))

# Compute all correlations between all variables considered for the model and EXT 
# as sanity check for the exclusion of variables from model
ext_vars <- indicators[indicators %in% names(df_corr)] # define externalizing variables in df_corr

# Compute correlations of each variable with EXT
# create  empty numeric vector to store correlations
correlations <- numeric(length(ext_vars)) 

# vector names to match variable names
names(correlations) <- ext_vars

# Loop through variables to get correltations with EXT
for (v in ext_vars) {
  correlations[v] <- cor( 
    df_corr[[v]], 
    df_corr$EXT_score, 
    use = "pairwise.complete.obs"
  )
}

# Put the results into a tibble for one row per variable
ext_vars_corr <- tibble(
  variable = names(correlations),
  correlation = correlations
) %>%
  arrange(desc(abs(correlation))) # order by corr

print(ext_vars_corr, n= 28)

# save as csv
write.csv(ext_vars_corr, file = file.path(path, "04_twinlife_ext_variables_correlations.csv"), row.names = FALSE)

# Create df for factor scores including demographic variables to export
df_fs <- as.data.frame(fs)

# Get row indices to map demographic data to factor scores for export
used_indices <- lavInspect(fit_ext_final, "case.idx")

# Map pid, fid, zygosity, sex and age in from data 
df_fs$pid <- df_model$pid[used_indices]
df_fs$fid <- df_model$fid[used_indices]
df_fs$zygosity <- df_model$zygosity[used_indices]
df_fs$sex <- df_model$sex[used_indices]
df_fs$age_yrs <- df_model$age_yrs_pq[used_indices]

# Check df
head(df_fs) 

# Count non-missing EXT scores
sum(!is.na(df_fs$EXT))  # 6851

# Check latent EXT score correlations with first oorder factors
cor_EXT <- cor(df_fs[, c("EXT", "hyp", "att", "srg", "con")], 
    use = "pairwise.complete.obs") # use pairwise to ignore NAs

cor_EXT

# Create heatmap EXT correlations using corrplot package
png(paste0(path, "/figures/", "EXT_Correlation_Heatmap.png"), width = 2000, height = 2000, res = 300) # prepare to save

corrplot(cor_EXT, method = "color", type = "upper", addCoef.col = "white", # plot corr 
         tl.col = "black", tl.srt = 45,                                  # set colors
         col = colorRampPalette(c("#8C8418", "#FFFFFF", "#5E1AF4"))(200),
         diag = FALSE, main = "Correlations Between EXT and Indicators") # title
dev.off()

# ============================================================
# 5. FACTOR COMPOSITION
# ============================================================

# Extract loadings for higher order EXT factor
std_res <- standardizedSolution(fit_ext_final, ci = TRUE) # standardized, with confidence intervals

# Get loadings for higher order EXT factor to plot
plot_data <- std_res[std_res$lhs == "EXT" & std_res$op == "=~", ]

# Create path diagram to show weight and structure of the model 
png(paste0(path, "/figures/", "EXT_path_diagram.png"), width = 2000, height = 2000, res = 300) # to save

# barplot
b_plot <- barplot(plot_data$est.std, # plot loadings
                  names.arg = plot_data$rhs,  # set x labels
                  col = "#5475E7",      
                  border = "white",       
                  main = "Composition of the Externalizing (EXT) Factor", # title
                  ylab = "Standardized loading (beta)", # y label
                  ylim = c(0, 1.1),       # leave room for error bars
                  cex.names = 1.1,        # make font bit larger
                  font.main = 2,
                  mar = c(5, 6, 4, 2))  # change margins for readability


# Add error bars as precision markers (95% CI)
arrows(b_plot, plot_data$ci.lower, b_plot, plot_data$ci.upper, # using 'lower' and 'upper' columns from standardizedSolution
       angle = 90, code = 3, length = 0.05, lwd = 1.5)

# Add values on top of bars 
text(x = b_plot, y = plot_data$est.std + 0.15,  # + 0.12 to avoid overlap error bars
     labels = round(plot_data$est.std, 2), 
     cex = 0.9, font = 2)

dev.off()

# Create another tree path diagram for whole factor composition 
png(paste0(path, "/figures/", "EXT_path_diagram2.png"), width = 2500, height = 2000, res = 300) # to save

semPaths(
  fit_ext_final,
  what = "std",          # standardized estimates
  whatLabels = "std",    # show loadings
  layout = "tree",       
  style = "lisrel",      # turned out better than others
  edge.color = "#4D4D4D",
  edge.width = 0.2, # thinner arrows
  
  color = list(
    man = "white",
    lat = c("white", "white", "white", "white", "#5E1AF4") # color EXT
  ),
  
  residuals = TRUE,
  intercepts = FALSE, # otherwise showed intercepts for observed variables also, messy
  edge.label.cex = 0.8, # smaller labels
  sizeMan = 6,
  sizeLat = 10,
  nCharNodes = 0,
  layoutScale = c(1.4, 1.2), # more space for labels
  optimizeLatRes = TRUE # push residuals from edges
)

# Add note on latent variance for EXT fixed to 1 
text(
  x = 0, # centered
  y = 1.5, # above EXT
  labels = "Var = 1"
)

dev.off()

# ============================================================
# 6. VARIANCE EXPLAINED + PCA
# ============================================================

# define variable with first order factor names
first_order_factors = c("hyp", "att", "srg", "con")

# 1. Get variance explained by EXT to check if there's a shared externalizing dimension
r2_vals <- inspect(fit_ext_final, "r2")
r2_first_order <- r2_vals[first_order_factors]

# How much higher-order EXT explains in each first-order factor
r2_first_order

# Variance EXT explains in first-order factors
mean(r2_first_order, na.rm = TRUE)


# 2. Inspect higher order loadings to check if EXT captures what first-order factors share
parameterEstimates(fit_ext_final, standardized = TRUE) %>%
  dplyr::filter(lhs == "EXT", op == "=~") %>%
  dplyr::select(lhs, rhs, est, std.all)

# Compare EXT to a composite of first-order factors
fs_compare <- df_fs %>%
  select(EXT, hyp, att, srg, con) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(c(hyp, att, srg, con), scale, .names = "{.col}_z")) %>%
  mutate(first_order_sum = hyp_z + att_z + srg_z + con_z )

# Overlap between EXT and composite
cor(fs_compare$EXT, fs_compare$first_order_sum, use = "pairwise.complete.obs") # correlation between EXT and standardized sum

# R^2 from regressing EXT on sum
summary(lm(EXT ~ first_order_sum, data = fs_compare))$r.squared 
  # higher-order score is almost identical to a standardized sum of first-order factor scores
  # supports that a general factor captures shared ranking, and some domains keep substantial unique variance


# 3. PCA as supporting evidence
# Create scree plot
png(paste0(path, "/figures/", "EXT_total_variance_explained.png"), width = 2000, height = 2000, res = 300)

# Define factor scores for pca
vars_pca <- df_fs[, first_order_factors]  

# Run diagnostic PCA on the complete rows
pca_diag <- prcomp(na.omit(vars_pca), scale. = TRUE)

# Calculate variance explained by first PC
var_explained <- (pca_diag$sdev^2 / sum(pca_diag$sdev^2))[1] * 100
var_explained


# Make scree Plot 
plot(pca_diag$sdev^2, type = "b", pch = 19, col = "#5E1AF4",
     main = "Scree Plot",
     xlab = "Component Number",
     ylab = "Variance Explained (Eigenvalue)",
     ylim = c(0, max(pca_diag$sdev^2) + 0.5),
     xaxt = "n")  # not draw default x acis

axis(side = 1, at = 1:length(pca_diag$sdev)) # define x axis with length of components
abline(h = 1, lty = 2, col = "#FF6230") # line at 1 kaiser criterion

dev.off()


# Variance explained by other PCs in percent
eig <- pca_diag$sdev^2 # eigenvalue
vars_explained <- eig / sum(eig) # variance proportional

pc2_var <- vars_explained[2] * 100
pc3_var <- vars_explained[3] * 100
pc4_var <- vars_explained[4] * 100
pc234_var <- sum(vars_explained[-1]) * 100 # sum of variance explained by all except first PC

# Print for reporting
cat("PC2 explains:", round(pc2_var, 2), "%\n")
cat("PC3 explains:", round(pc3_var, 2), "%\n")
cat("PC4 explains:", round(pc4_var, 2), "%\n")
cat("Sum PC2-4:" , round(pc234_var, 2), "%\n") # PC1 explains waay more than sum of other PCs

# inspect PCs and first order factors
round(pca_diag$rotation, 2)

# Add first PC per pid to df_fs
# Variables used for PCA
pc_vars <- first_order_factors

# Identify rows with complete data on PCA variables
pc_complete <- complete.cases(df_fs[, pc_vars])

# Run PCA on complete cases only
pca_diag <- prcomp(
  df_fs[pc_complete, pc_vars],
  center = TRUE,
  scale. = TRUE
)

# create empty PC1 numeric column first
df_fs$PC1 <- NA_real_

# Insert PC1 back into matching rows
df_fs$PC1[pc_complete] <- pca_diag$x[, 1]

# Sanity check
summary(df_fs$PC1)
cor(df_fs$EXT, df_fs$PC1, use = "pairwise.complete.obs") # highly correlated
                                                        # but also higher PC1 currently indicates lower externalizing !

# Reverse PC1 so higher PC1 indicates higher externalizing
df_fs$PC1 <- -df_fs$PC1

# Check reverse coding
cor(df_fs$EXT, df_fs$PC1, use = "pairwise.complete.obs") # now higher PC1 = higher externalizing

# Standardize PC1 for downstream analyses
df_fs$PC1 <- as.numeric(scale(df_fs$PC1))

# Check standardization
mean(df_fs$PC1, na.rm = TRUE)  # close to 0, good
sd(df_fs$PC1, na.rm = TRUE)    # 1, perfect

# PC1 summary statistics
summary(df_fs$PC1)

# Standardize EXT and first order factor scores
df_fs <- df_fs %>%
  mutate(
    across(
      all_of(c(first_order_factors, "EXT")),
      ~ as.numeric(scale(.x))
    )
  )


# ============================================================
# 7. SAVE DATA
# ============================================================

# Remove not used variables to save final model data
df_model <- df_model %>%
  select(
    all_of(demo_vars),
    all_of(indicators)
  )

# Save model data
write.csv(df_model, file.path(path, "04_twinlife_externalizing_model_data.csv"), row.names = FALSE)
save(df_model, file = file.path(path, "04_twinlife_externalizing_model_data.rda"))

# Save factor scores and PC1 including sex, zygosity, age in years 
write.csv(df_fs, file.path(path, "04_twinlife_externalizing_cfa.csv"), row.names = FALSE)
save(df_fs, file = file.path(path, "04_twinlife_externalizing_cfa.rda"))
