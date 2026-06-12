########################################################################
#### Externalizing — Factor Analysis **Model 3** in TwinLife        #### 
#### Sensitivity model with biggest Modification Index              ####
#### 01.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model3 <- '
  # ADHD related indicators
  hyp =~ adhd_restless + adhd_move 
  att =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation
  srg =~ srg0100 + srg0200 + srg0400 + srg0500 + srg0600
  
  # Conduct problems: disruptive behavior and delinquency
  con =~ agg + devi_anger + del_lie + del_steal
  
  # Substance use
  sub =~ smok + alc_freq_exc + alc_n_beer_exc + alc_n_hp_exc

  # Higher-order factor
  EXT =~ hyp + att + srg + con + sub    
         
  # Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  srg0100 ~~ srg0200

  # Sensitivity model including srg with alc_freq_excand smok residual covariance
  srg0400 ~~ alc_freq_exc
  srg0400 ~~ smok
'