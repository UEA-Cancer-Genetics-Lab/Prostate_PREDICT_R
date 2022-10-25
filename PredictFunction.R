# Implementation of the PREDICT Prostate clinical prognostic model in R
# Adapted from the STRATA implementation from Thurtle DR, Greenberg DC, Lee LS, Huang HH, Pharoah PD, Gnanapragasam VJ. Individual prognosis at diagnosis in nonmetastatic prostate cancer: Development and external validation of the PREDICT Prostate multivariable model. Johnstone RW, editor. PLoS Medicine [Internet]. 2019 Mar;16(3):e1002758. Available from: http://dx.plos.org/10.1371/journal.pmed.1002758

library(tidyverse)
library(reshape2)

# NPCM PCSM is the probability for nonPCa and PCa death respectively
# i is years we want to predict

predict_function <- function(i,
                             age,
                             gradegroup,
                             psa,
                             t_stage,
                             charlson_comorbidity,
                             primaryRx,
                             biopsy50){
  
  time_y <- seq(1,i) 
  
  dat_1 <- seq(1,i) %>%
    as.data.frame() %>%
    dplyr::rename(time = ".") %>%
    mutate(age = age) %>%
    mutate(
      age_time = age + time,
      gradegroup= gradegroup,
      psa = psa, # ng/ml
      t_stage =t_stage,
      charlson_comorbidity = charlson_comorbidity, # charlson comorbidity score (simplified to 0 or 1 (1 or greater)
      primaryRx = primaryRx, # 0=Conservative Mx 1='radical treatment' 3=ADT alone
      biopsy50 = biopsy50) # 0=Unknown/not included 1=<50% cores involved 2=>=50% cores involved (Enter 0 if using only the baseline model)
  
  
  # calculate the PCSM prognostic index (pi)
  dat_2 <- dat_1 %>% 
    mutate(piPCSM = 0.0026005*((age/10)^3-341.155151) +  # calculate the PCSM prognostic index (pi)
             0.185959*(log((psa+1)/100)+1.636423432)  + 
             .1614922*(t_stage==2) + .39767881*(t_stage==3) + 
             .6330977*(t_stage==4) + .2791641*(gradegroup==2) + 
             .5464889*(gradegroup==3) + .7411321*(gradegroup==4) + 
             1.367963*(gradegroup==5) + -.6837094*(primaryRx==1) + 
             .9084921*(primaryRx==3) -0.617722958*(biopsy50==1) + 
             0.579225231*(biopsy50==2)) %>%
    mutate(piNPCM = 0.1226666*(age-69.87427439) + 0.6382002*(charlson_comorbidity==1)) %>% #calculate the NPCM progostic index (pi)
    rowwise() %>%
    mutate(time = time * 365) %>% # change to days
    mutate(PCSMatT = 1 - exp(-exp(piPCSM)*exp(-16.40532 + 1.653947*(log(time)) + 1.89e-12*(time^3)))) %>% #PCS mortality, then per year, then convert to survival per year...
    mutate(NPCMatT = 1 - exp(-exp(piNPCM)*exp(-12.4841 + 1.32274*(log(time)) + 2.90e-12*(time^3))))
  
  # PCS mortality, then per year, then convert to survival per year...
  PCSMatT <- dat_2$PCSMatT
  PCSMatT_sub <- c(0,dat_2$PCSMatT[1:length(time_y)-1])
  PCSM_mortrate_year <- PCSMatT - PCSMatT_sub
  PCSsurvival_year <- 1 - PCSM_mortrate_year
  
  dat_3 <- cbind(dat_2,PCSM_mortrate_year,PCSsurvival_year)
  
  
  # Do the same for NPCM...
  NPCMatT <- dat_2$NPCMatT
  NPCMatT_sub <- c(0,dat_2$NPCMatT[1:length(time_y)-1])
  NPCMatT_mortrate_year <- NPCMatT - NPCMatT_sub
  NPCsurvival_year <- 1 - NPCMatT_mortrate_year
  
  dat_4 <- cbind(dat_3,NPCMatT_mortrate_year,NPCsurvival_year) %>%
    mutate(PCSsurvival = 1 - PCSMatT) %>%
    mutate(NPCsurvival = 1 - NPCMatT) %>%
    mutate(allcauseM = 1 - PCSsurvival*NPCsurvival)
  
  # 'generate all cause mortality'
  allcauseM <- dat_4$allcauseM
  allcauseM_sub <- c(0,dat_4$allcauseM[1:length(time_y)-1])
  allcauseM_inyear <- allcauseM - allcauseM_sub
  
  dat_5 <- cbind(dat_4,allcauseM_inyear) %>%
    mutate(proportionPC_cum = PCSMatT / (PCSMatT + NPCMatT)) %>% #calculate proportion of all cause mortality from cause-specific death
    mutate(proportionPC = (PCSM_mortrate_year)/(NPCMatT_mortrate_year + PCSM_mortrate_year)) %>%
    mutate(propn_NPC = 1-proportionPC) %>% #PC mortality as competing risk
    mutate(pred_PC_year = proportionPC*allcauseM_inyear)
  
  pred_PC_year = dat_5$pred_PC_year
  pred_PC_cum <- cumsum(pred_PC_year)
  
  dat_6 <- cbind(dat_5,pred_PC_cum) %>%
    mutate(pred_NPC_year = propn_NPC*allcauseM_inyear)
  
  pred_NPC_cum <- dat_6$pred_NPC_year
  pred_NPC_cum <- cumsum(pred_NPC_cum)
  
  dat_7 <- cbind(dat_6,pred_NPC_cum) %>%
    mutate(overallsurvival = 100) %>%
    mutate(NPCM = pred_NPC_cum*100 + pred_PC_cum*100) %>%
    mutate(PCSM = pred_PC_cum*100)
  
  return(dat_7)
}
