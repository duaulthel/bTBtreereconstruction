rm(list=ls())

#--------------------------------
#packages
library(tidyverse) #version 1.3.0

#--------------------------------
#Description:
# Check for the presence of super-spreaders
# Definition of super-spreaders: less than 10% of infected hosts 
# are responsible for over 80% of transmission events

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes

#--------------------------------

for (j in 1:30){ #j is the tree number
  
  sim <- paste0(samp,"_",j)
  
  #initialize output
  pres_spp <- NULL
  
  for (i in 1:nb_scheme){#sampling schemes
    
    #results from seqTrack
    get_res <- read_csv(paste0("./seqTrack_",samp,"/seqTrack_tree",sim,"_",i,".csv")) 

    #remove index case (infector == NA) in order to consider transmission events
    data_sic <- get_res %>% filter(!is.na(infector))
    
    #Calculate number of secondary cases per infector
    nb_sc <- data_sic %>% dplyr::count(ances)
    
    ##Presence of super-spreaders?
    sup <- tail(nb_sc %>% arrange(n), ifelse(floor(nrow(data_sic)*0.1) >0, floor(nrow(data_sic)*0.1), 1))
    #select the 10% of infected hosts that have the most secondary cases
    
    #Create output for this tree
    correct <- tibble("sim"=sim,
                      "scenario"= i, 
                      "ssp"=ifelse(sum(sup$n)>(0.8*nrow(data_sic)), 1, 0)) 
    #ssp notes the presence or not of superspreaders in the reconstructed tree
    
    pres_spp <- rbind(pres_spp, correct)
  }

  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("seqTrack_",samp,"_ssp.csv"))
    pres_spp <- rbind(prev, pres_spp)
    write_csv(pres_spp, paste0("seqTrack_",samp,"_ssp.csv"))
  }
  else{ #if this is the first tree, write new file
    write_csv(pres_spp, paste0("seqTrack_",samp,"_ssp.csv"))
  }
  
}
