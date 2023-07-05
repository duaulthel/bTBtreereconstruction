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
samp <- "B1"

nb_scheme <- ifelse(samp=="B1", 6, 1)
#------------------------------------

for (j in 1:30){ #j is the tree number
  
  sim <- paste0(samp,"_",j)
  
  #find index case of the reference tree
  real_index <- ifelse(substr(samp, 2,2)=="1", "cattle", "badger")
  
  #initialize output
  index <- NULL

  for (i in 1:nb_scheme){ #sampling schemes

    #sampled sequences that contains host-species information
    seq <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",sim,"_",i,".csv"))
    seq$id <- 1:nrow(seq) #add id used in outbreaker2

    #results from outbreaker2
    tree <- read_csv(paste0("./outbreaker_",samp,"/res_tree",sim,"_",i,".csv")) 

    #merge tree with sequences with host-species information
    #remove nucleotide column in seq (number 4)
    data <- merge(seq[,-4], tree, by.x="id", by.y="to")
    
    #find the name of infectors and their host-species
    data <- data %>% mutate(infector=infected[from],
                            sp_infector=sp[from])
    
     #Find index case's host-species
    data_ic <- data %>% filter(!is.na(support)) %>%
      filter(time==min(time)) %>% 
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
    prev <- read_csv(paste0("outbreaker_",samp,"_index.csv"))
    index <- rbind(prev, index)
    write_csv(index, paste0("outbreaker_",samp,"_index.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(index, paste0("outbreaker_",samp,"_index.csv"))
  }
  
}
