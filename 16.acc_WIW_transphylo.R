rm(list=ls())
#------------------------------------
#package
library(tidyverse)

#------------------------------------
##Description:
##Estimate number of correctly reconstructed transmission
##events and number of incorrectly reconstructed events.

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#Transmission scenario considered
samp <- "B1"

#------------------------------------
non_se <- c("A1", "B1", "B2", "S1") 
#tree_name of reference tree same as reconstructed tree name

#For the reconstructed trees with a higher mutation rate
#Find name of the reference tree
ref_tree <- case_when(samp %in% non_se ~ samp,
                      samp=="S4" ~ "B1")

nb_scheme <- ifelse(samp=="B1", 6, 1)
#------------------------------------

#Not all trees converged, list of those that did:
# samp=="A1" ~ c(1,6,7,8,10,13,16,17,20,21,22,24,25,26,27,28,29)
# samp=="B1" ~ c(1:4,8,10:14,17,18,20:23,26:30)
# samp=="B2" ~ c(1:6,9,15:18,20,22:26,29)
# samp=="S1" ~ c(1:2,4,6:15,17:26,28:30)
# samp=="S4" ~ c(1:3,6,8:9,11:15,17:18,20:22,24:27,29)

#----------------------------------------------------------------------------

for (j in c(1:2,4:6,8:9,11:16,18:19,22:29)){
  
  sim <- paste0(samp,"_",j) 
  
  #reference tree with details from Ttree_biased_trans
  Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_det_",ref_tree,"_",j,".csv")) #reference tree

  #remove index case (infector == NA) in order to consider transmission events 
  Ttree <- Ttree %>% filter(!is.na(infector))
  
  #initialize accuracy output
  acc <- NULL
  
  for (i in 1:nb_scheme){ #sampling schemes

    #sampled sequences with host-species information
    seq <- read_csv(paste0("./seq_",substr(samp, 1, 1),"/seq_",sim,"_",i,".csv"))
    tstart <- as.Date("01/01/07", "%d/%m/%y")
    seq$tremoved <- tstart %m+% months(seq$tremoved)
    seq <- seq %>% 
      mutate(infected_id = paste(infected,sp,tremoved, sep="_")) 
    #in order to have the same name as in beast
    
    #results from matrix TransPhylo
    tree <- read_csv(paste0("./TransPhylo_",samp,"/WIWtree_",sim,"_",i,".csv")) 
    
    #add host-species information for infected
    data <- full_join(seq[,c(1,3,5)], tree, by="infected_id")
    
    #Rename infected into infector in seq
    seq <- seq %>% rename(infector=infected, sp_infector=sp, infector_id=infected_id)
    
    #add host-species information for infectors
    data <- right_join(seq[,c(1,3,5)], data, by="infector_id")
    
    #Create new columns
    data_sic <- data %>% filter(prob > 0) %>% #remove those whose infector is unknown
      mutate(trans = paste(infector, infected, sep="->"),
             transmission = paste(sp_infector, sp, sep="->"))
    
    #Proportion of correctly estimated transmission events
    correct <- data_sic %>% filter(trans %in% Ttree[[paste0("trans_",i)]])
    
    #Create tibble with results for accuracy
    out <- tibble("sim"=sim,
                  "scenario"= i,
                  "nb_acc"=nrow(correct),
                  #number of 
                  #correctly estimated transmission events
                  "nb_non_acc"=nrow(data_sic)-nrow(correct))
    
    acc <- rbind(acc, out)
  }
  
  #------------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("transphylo_WIW_",samp,"_new_acc.csv"))
    acc <- rbind(prev, acc)
    write_csv(acc, paste0("transphylo_WIW_",samp,"_new_acc.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(acc, paste0("transphylo_WIW_",samp,"_new_acc.csv"))
  }
  
}