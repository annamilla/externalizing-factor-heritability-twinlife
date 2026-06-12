########################################################################
#### Externalizing — Factor Analysis **Model 1** in TwinLife        #### 
#### Initial model including all indicators with good data coverage ####
#### 01.04.2026 - 12.06.2026                                        ####
#### main Author: Anna Miller                                       #### 
########################################################################

model0 <- '
  # ADHD related indicators
  hyp =~ adhd_restless + adhd_move 
  att =~ adhd_act + adhd_tasks + adhd_unfoc
  
  # Self-regulation
  srg =~ srg0100 + srg0200 + srg0300 + srg0400 + srg0500 + srg0600
  
  # Disruptive behavior
  dis =~ agg + devi_anger
  
  # Delinquency
  del =~ del_compl + del_lie + del_steal 
  
  # Substance use
  sub =~ smok + alc_freq_exc + alc_init + alc_init_exc + alc_n_beer_exc +
        alc_n_hp_exc + alc_n_wine_exc 

  # Higher-order externalizing factor
  EXT =~ hyp + att + srg + dis + del + sub    
'