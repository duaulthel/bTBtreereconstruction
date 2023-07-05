rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(lubridate)

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

#If there is a sampling scheme different from A, we need the tree from induced_subtree
sup <- ifelse(substr(samp,1,1)=="S", "det", "det_sup")

#--------------------------------
#Not all trees converged, list of those that did:
# samp=="A1" ~ c(1,6,7,8,10,13,16,17,20,21,22,24,25,26,27,28,29)
# samp=="B1" ~ c(1:4,8,10:14,17,18,20:23,26:30)
# samp=="B2" ~ c(1:6,9,15:18,20,22:26,29)
# samp=="S1" ~ c(1:2,4,6:15,17:26,28:30)
# samp=="S4" ~ c(1:3,6,8:9,11:15,17:18,20:22,24:27,29)

#----------------------------------------------------------------------------

for (j in c(1:2,4:6,8:9,11:16,18:19,22:29)){#j is the tree considered
  
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

    #results from medoid TransPhylo
    tree <- read_csv(paste0("./TransPhylo_",samp,"/ttree_",sim,"_",i,".csv"))
    
    #Create tibble with results from size estimations
    nb_comp <- tibble(nb_ref=nrow(Ttree), #size of the induced subtree
             nb_sim=nrow(tree), #size of the reconstructed tree
             tree=sim, 
             scenario=i)
    
    nb <- rbind(nb, nb_comp)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("transphylo_",samp,"_tree_size.csv"))
    nb <- rbind(prev, nb)
    write_csv(nb, paste0("transphylo_",samp,"_tree_size.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(nb, paste0("transphylo_",samp,"_tree_size.csv"))
  }
  
}
