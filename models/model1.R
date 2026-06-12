########################################################################
#### Externalizing — Factor Analysis **Model 1** in TwinLife        #### 
#### Model including substance use as first order factor            ####
#### 01.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model1 <- '
  # ADHD related indicators
  hyp =~ adhd_restless + adhd_move 
  att =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation 
  srg =~ srg0100 + srg0200 + srg0400 + srg0500 + srg0600    # removed srg0300 with low coverage
  
  # Conduct problems: disruptive behavior and delinquency 
  con =~ agg + devi_anger + del_lie + del_steal       # removed del_compl with low coverage
  
  # Substance use 
  sub =~ smok + alc_freq_exc + alc_n_beer_exc + alc_n_hp_exc    # removed alc_init_exc, alc_n_wine_exc with low coverage

  # Higher-order factor
  EXT =~ hyp + att + srg + con + sub    
         
  #  Add residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  adhd_restless ~~ adhd_move
  srg0100 ~~ srg0200              
'