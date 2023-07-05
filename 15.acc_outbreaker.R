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
for (j in 1:30){#all trees in transmission scenario
  
  sim <- paste0(samp,"_",j) 
  
  #reference tree with details from Ttree_biased_trans
  Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_det_",ref_tree,"_",j,".csv")) 
  
  #remove index case (infector == NA) in order to only consider transmission events 
  Ttree <- Ttree %>% filter(!is.na(infector))
  
  #initialize accuracy output
  acc <- NULL
  
  for (i in 1:nb_scheme){ #sampling schemes

    #sampled sequences with host-species information
    seq <- read_csv(paste0("./seq_",substr(samp, 1, 1),"/seq_",sim,"_",i,".csv")) 
    seq$id <- 1:nrow(seq)#add id used in outbreaker2

    #results from outbreaker2
    tree <- read_csv(paste0("./outbreaker_",samp,"/res_tree",sim,"_",i,".csv")) 
    
    #merge tree with sequences with host-species information
    data <- merge(seq[,-4], tree, by.x="id", by.y="to")
    
    #find the name of infectors and their host-species
    data <- data %>% mutate(infector=infected[from],
                            sp_infector=sp[from])
    
    #remove index case (generations == NA) in order to consider direct transmission events
    data_sic <- data %>%
      filter(!is.na(generations)) %>%
      mutate(trans = paste(infector, infected, sep="->"),
             transmission = paste(sp_infector, sp, sep="->"))

    #Correctly estimated transmission events
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
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("outbreaker_",samp,"_acc.csv"))
    acc <- rbind(prev, acc)
    write_csv(acc, paste0("outbreaker_",samp,"_acc.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(acc, paste0("outbreaker_",samp,"_acc.csv"))
  }
  
}
