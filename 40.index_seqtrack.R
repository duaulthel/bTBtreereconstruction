rm(list=ls())

#--------------------------------
#packages
library(tidyverse)

#--------------------------------
#Description:
## Identify whether the correct index case's host-species 
## was reconstructed for each tree

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#Transmission scenario considered
samp <- "B1"

nb_scheme <- ifelse(samp=="B1", 6, 1)
#------------------------------------

for (j in 1:30){#j is the tree number
  
  sim <- paste0(samp,"_",j)
  
  #find index case of the reference tree
  real_index <- ifelse(substr(samp, 2,2)=="1", "cattle", "badger")
  
  #initialize output
  index <- NULL
  
  for (i in 1:nb_scheme){ #sampling schemes
    
    #results from seqTrack
    get_res <- read_csv(paste0("./seqTrack_",samp,"/seqTrack_tree",sim,"_",i,".csv")) 

    #Find index case's host-species
    data_ic <- get_res %>% 
      filter(is.na(infector))%>% 
      filter(date==min(date)) %>%  
      group_by(sp) %>%
      summarise(nb=n_distinct(infected)) 
    #count number of index case's per host-species
    
    #is there a host-species more frequently reconstructed as index
    maj_index <- data_ic %>% slice_max(nb)
    
    #Identify most frequently reconstructed host-species
    #If there isn't one, index_sim equals 0
    index_sim <-  ifelse(nrow(maj_index==1), maj_index$sp, 0)

    #Create output for this tree
    correct <- tibble("sim"=sim,
                      "scenario"= i, 
                      "index"=ifelse(index_sim==real_index, 1, 0)) #index case's host-species
    
    index <- rbind(index, correct)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("seqTrack_",samp,"_index.csv"))
    index <- rbind(prev, index)
    write_csv(index , paste0("seqTrack_",samp,"_index.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(index , paste0("seqTrack_",samp,"_index.csv"))
  }
  
}