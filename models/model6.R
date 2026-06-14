########################################################################
#### Externalizing — Factor Analysis **Model 6** in TwinLife        #### 
#### Sensitivity check one order model                              ####
#### 30.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model6 <- '
  # Externalizing Factor for all indicators
  EXT =~ adhd_restless + adhd_move + 
         adhd_act + adhd_tasks + adhd_unfoc +
         srg0100 + srg0200 + srg0400 + srg0500 + srg0600 +
         agg + devi_anger + del_lie + del_steal
  
         
  # Residual covariances to account for method effects
  adhd_unfoc ~~ adhd_tasks
  srg0100 ~~ srg0200              
'


