########################################################################
#### Externalizing — Factor Analysis **Model 2** in TwinLife        #### 
#### Model without residual cov restless/move                       ####
#### 01.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model2 <- '
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

  #  Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks        # removed residual cov restless/move for two-indicator factor hyp
  srg0100 ~~ srg0200              
'
