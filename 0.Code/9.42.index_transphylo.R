rm(list=ls())

#--------------------------------
#packages
library(tidyverse) #version 1.3.0
library(lubridate) #version 1.7.9.2

#--------------------------------
#Description:
## Identify whether the correct index case's host-species 
## was reconstructed for each tree

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1" #B1 (reference scenario) or B2 (badger index)

nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes
#--------------------------------

#Not all trees converged, list of those that did:
# samp=="B1" ~ c(1:4,8,10:14,17,18,20:23,26:30)
# samp=="B2" ~ c(1:6,9,15:18,20,22:26,29)

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

for (j in c(1:4,8,10:14,17,18,20:23,26:30)){ #change to list of trees that converged
  
  sim <- paste0(samp,"_",j)
  
  #find index case of the reference tree
  real_index <- ifelse(substr(samp, 2,2)=="1", "cattle", "badger")
  
  #initialize output
  index <- NULL 
  
  for (i in 1:nb_scheme){ #sampling schemes
  
    #sampled sequences with host-species information
    seq <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",sim,"_",i,".csv"))
    tstart <- as.Date("01/01/07", "%d/%m/%y") #start of the outbreak
    seq$tremoved <- tstart %m+% months(seq$tremoved)
    seq <- seq %>% 
      mutate(infected_id = paste(infected,sp,tremoved, sep="_")) 
    #create infected_id in order to have the same name as in beast
    
  
    #results from TransPhylo_medoid
    tree <- read_csv(paste0("./TransPhylo_",samp,"/ttree_",sim,"_",i,".csv")) 
    
    #find the name of infectors and their host-species
    tree <- tree %>% mutate(infector_id=ifelse(infector_id!=0,infected_id[infector_id], NA))
    
    #merge tree and seq in order to get host-species information
    data <- full_join(seq[,c(1,3,5)], tree, by="infected_id")
    
    #Rename infected into infector in seq
    seq <- seq %>% rename(infector=infected,
                          sp_infector=sp,
                          infector_id=infected_id)
    
    #merge tree and seq in order to get host-species information
    data <- right_join(seq[,c(1,3,5)], data, by="infector_id")
    
    #Find index case's host-species
    index_sim <- data[which.min(data$tinfection),]$sp
    index_sim <- ifelse(is.na(index_sim), "NA", index_sim)

    #Create output for this tree
    correct <- tibble("sim"=sim,
                      "scenario"= i, 
                      "index"=ifelse(index_sim==real_index, 1, 0)) #index case's host-species
    
    index <- rbind(index, correct)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("transphylo_",samp,"_index.csv"))
    index <- rbind(prev, index)
    write_csv(index, paste0("transphylo_",samp,"_index.csv"))
  }
  else{ #if this is the first tree, write new file
    write_csv(index, paste0("transphylo_",samp,"_index.csv"))
  }
  
}
