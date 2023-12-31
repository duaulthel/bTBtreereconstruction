rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(reshape2) #version 0.8.8

#------------------------------------
##Description:
##Transform matrix from TransPhylo into readable dataframe
##with the most likely infector 

#------------------------------------
#------------------------------------
#Transmission scenario
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes

#------------------------------------
#Function that transforms matrix into a tree with the most likely infector
mat_to_ttree <- function(WIW){
  rownames(WIW) <- colnames(WIW) 
  tree <- setNames(melt(as.matrix(WIW)), c('infector_id', 'infected_id', 'prob'))
  tree <- tree %>% group_by(infected_id) %>% slice(which.max(prob))
  return(tree)
}

#------------------------------------
#apply function mat_to_ttree to TransPhylo results

for (j in 1){ #trees that converged
  for (i in 1:nb_scheme){ #sampling schemes
    sim <- paste0(samp,"_",j)
    
    #results from TransPhylo
    WIW <- read_csv(paste0("./TransPhylo_",samp,"/WIW_",sim,"_",i,".csv")) 
    
    #apply function
    tree <- mat_to_ttree(WIW)
    
    #Write output
    if (!is.null(tree)){
      write_csv(tree, paste0("WIWtree_",sim,"_",i,".csv"))
    }
  }
}



