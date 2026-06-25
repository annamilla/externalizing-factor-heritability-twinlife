########################################################################
#### Externalizing — Factor Analysis **Model Final** in TwinLife    #### 
#### Final Model without substance use                              ####
#### 30.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model_final <- '
  # ADHD related indicators
  HYP =~ adhd_restless + adhd_move
  ATT =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation
  SCT =~ srg0100 + srg0200 + srg0400 + srg0500 + srg0600
  
  # Conduct problems: disruptive behavior and delinquency
  CON =~ agg + devi_anger + del_lie + del_steal

  # Higher-order factor without substance use
  EXT =~ HYP + ATT + SCT + CON  
         
  # Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  srg0100 ~~ srg0200              
'


