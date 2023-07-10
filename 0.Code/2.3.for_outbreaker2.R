rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(epitrix) #version 0.2.2
library(lubridate) #version 1.7.9.2
library(fitdistrplus) #version 1.1-6

#------------------------------------
##Description:
##Create distributions for generation times and sampling intervals 
##needed to implement outbreaker2 

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "A1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)


for (j in 1:30){ #all trees in the transmission scenario
  
  Ttree <- read_csv(paste0("./Ttrees_",substr(samp,1,1),"/Ttree_",samp,"_",j,".csv")) #reference tree
  
  #------------------------------------
  #Generation time distribution
  
  #Get all generation times from reference tree
  out <- NULL
  for(l in Ttree$infected[-1]){
    sub_tree <- Ttree %>% filter(infected == l) 
    gen_inf <- sub_tree$infector  
    gen_time <- sub_tree$tinfection-Ttree[Ttree$infected == gen_inf,]$tinfection
    out <- c(out, gen_time)
  }

  #For outbreaker: create distribution
  distr_gen <- as.vector(table(out))
  distr_gen <- distr_gen*1/(length(out))
  
  #Write output for generation times
  write.csv(distr_gen, paste0("generation_time_",samp,"_",j,".csv")) 
  
  #------------------------------------
  #Sampling time distribution
  
  #Get all sampling intervals from reference tree
  out <- NULL
  for(k in Ttree$infected){
    sub_tree <- Ttree %>% filter(infected == k) 
    samp_time <- sub_tree$tremoved-sub_tree$tinfection
    out <- c(out, samp_time)
  }

  #For outbreaker: create distribution
  distr_samp <- as.vector(table(out))
  distr_samp <- distr_samp*1/(length(out))
  
  #Write output for sampling intervals
  write.csv(distr_samp, paste0("sampling_time_",samp,"_",j,".csv")) 
}
