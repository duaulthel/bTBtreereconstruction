rm(list=ls())
#------------------------------------
#package
library(tidyverse)

#------------------------------------
##Description:
##Estimate size of the induced subtree
##Estimate size of the reconstructed tree

#------------------------------------
#------------------------------------
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

#If there is a sampling scheme different form A, we need the tree from induced_subtree
sup <- ifelse(substr(samp,1,1)=="S", "det", "det_sup")
#--------------------------------

for (j in 1:30){#j is the tree considered
  
  #initialize output
  nb <- NULL
  
  for (i in 1:nb_scheme){#sampling schemes
    
    sim <- paste0(samp,"_",j) 
    
    #reference tree
    Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_",sup,"_",ref_tree,"_",j,".csv")) 

    #if sampled scheme, take into account the induced subtree
    #> remove hosts not present in the subtree
    if(i>1){
      Ttree <- Ttree[(Ttree[[paste0("presence_",i)]])==1,]
    }
    
    
    #results from outbreaker2
    tree <- read_csv(paste0("./outbreaker_",samp,"/res_tree",sim,"_",i,".csv")) 
    
    #Remove index case to sum generations
    #generation = 1 no unsampled intermediary
    data_sic <- tree %>% filter(!is.na(generations))
    
    #Estimate size of reconstructed tree by accounting for unsampled hosts
    nb_hosts <- nrow(tree) + (sum(data_sic$generations)-nrow(data_sic))
    
    #Create tibble with results from size estimations
    nb_comp <- tibble(nb_ref=nrow(Ttree), #size of the induced subtree
                      nb_sim=nb_hosts, #size of the reconstructed tree
                      tree=sim, 
                      scenario=i)
    
    nb <- rbind(nb, nb_comp)
    
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("outbreaker_",samp,"_tree_size.csv"))
    nb <- rbind(prev, nb)
    write_csv(nb, paste0("outbreaker_",samp,"_tree_size.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(nb, paste0("outbreaker_",samp,"_tree_size.csv"))
  }
  
}
