########################################################################
#### Externalizing — Factor Analysis **Model 5** in TwinLife        #### 
#### test whether substance use is related to EXT                   ####
#### 30.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model5 <- '
  # ADHD related indicators
  hyp =~ adhd_restless + adhd_move
  att =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation
  srg =~ srg0100 + srg0200 + srg0400 + srg0500 + srg0600
  
  # Conduct problems: disruptive behavior and delinquency
  con =~ agg + devi_anger + del_lie + del_steal

  # Substance use as separate factor
  sub =~ smok + alc_freq_exc + alc_n_beer_exc + alc_n_hp_exc

  # Higher-order factor without substance use
  EXT =~ hyp + att + srg + con    

  # Association between EXT and substance use factors
  EXT ~~ sub
         
  # Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  srg0100 ~~ srg0200              
'

