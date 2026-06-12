########################################################################
#### Externalizing — Factor Analysis **Model 4** in TwinLife        #### 
#### Sensitivity model with composite instead of latent factors     ####
#### Test if the 2-indicator factor is driving higher-order model   ####
#### 01.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model4 <- '
  # ADHD related indicators without hyp
  att =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation
  srg =~ srg0100 + srg0200 + srg0400 + srg0500 + srg0600
  
  # Conduct problems: disruptive behavior and delinquency
  con =~ agg + devi_anger + del_lie + del_steal
  
  # Substance use
  sub =~ smok + alc_freq_exc + alc_n_beer_exc + alc_n_hp_exc

  # Higher-order factor including composite hyp score instead of latent hyp factor
  EXT =~ att + hyp_comp + srg + con + sub

  # Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  srg0100 ~~ srg0200
'


