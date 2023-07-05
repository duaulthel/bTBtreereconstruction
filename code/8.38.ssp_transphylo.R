rm(list=ls())

#--------------------------------
#packages
library(tidyverse)

#--------------------------------
#Description:
# Check for the presence of super-spreaders
# Definition of super-spreaders: less than 10% of infected hosts 
# are responsible for over 80% of transmission events

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "A1"

nb_scheme <- ifelse(samp=="B1", 6, 1)
#--------------------------------

#Not all trees converged, list of those that did:
# samp=="A1" ~ c(1,6,7,8,10,13,16,17,20,21,22,24,25,26,27,28,29)
# samp=="B1" ~ c(1:4,8,10:14,17,18,20:23,26:30)
# samp=="B2" ~ c(1:6,9,15:18,20,22:26,29)
# samp=="S1" ~ c(1:2,4,6:15,17:26,28:30))
# samp=="S4" ~ c(1:3,6,8:9,11:15,17:18,20:22,24:27,29)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

for (j in c(1:6,9,15:18,20,22:26,29)){ #j is the tree number
  
  sim <- paste0(samp,"_",j)
  
  #initialize output
  pres_spp <- NULL

  for (i in 1:nb_scheme){#sampling schemes
    
    #results from TransPhylo_medoid
    tree <- read_csv(paste0("./TransPhylo_",samp,"/ttree_",sim,"_",i,".csv")) 
    
    #Calculate number of secondary cases per infector
    nb_sc <- tree %>% dplyr::count(infector_id)
    
    ##Presence of super-spreaders?
    sup <- tail(nb_sc %>% arrange(n), ifelse(floor(nrow(tree)*0.1) >0, floor(nrow(tree)*0.1), 1))
    #select the 10% of infected hosts that have the most secondary cases
    
    #Create output for this tree
    correct <- tibble("sim"=sim,
                      "scenario"= i, 
                      "ssp"=ifelse(sum(sup$n)>(0.8*nrow(tree)-1), 1, 0)) 
    #ssp notes the presence or not of superspreaders in the reconstructed tree
    
    pres_spp <- rbind(pres_spp, correct)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("transphylo_",samp,"_ssp.csv"))
    pres_spp <- rbind(prev, pres_spp)
    write_csv(pres_spp, paste0("transphylo_",samp,"_ssp.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(pres_spp, paste0("transphylo_",samp,"_ssp.csv"))
  }
  
}
