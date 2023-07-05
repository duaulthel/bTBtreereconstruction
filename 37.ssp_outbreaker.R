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
samp <- "B1"

nb_scheme <- ifelse(samp=="B1", 6, 1)

#--------------------------------

for (j in 1:30){#j is the tree number
  
  sim <- paste0(samp,"_",j)
  
  #initialize output
  pres_spp <- NULL
  
  for (i in 1:nb_scheme){#sampling schemes

    #results from outbreaker2
    tree <- read_csv(paste0("./outbreaker_",samp,"/res_tree",sim,"_",i,".csv")) #results from outbreaker2
    
    #remove index case (generations == NA) in order to consider direct transmission events
    data_sic <- tree %>% filter(!is.na(generations))
    
    #Calculate number of secondary cases per infector
    nb_sc <- data_sic %>% dplyr::count(from)
    
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
    prev <- read_csv(paste0("outbreaker_",samp,"_ssp.csv"))
    pres_spp <- rbind(prev, pres_spp)
    write_csv(pres_spp, paste0("outbreaker_",samp,"_ssp.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(pres_spp, paste0("outbreaker_",samp,"_ssp.csv"))
  }
  
}
