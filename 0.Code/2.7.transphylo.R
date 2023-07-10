rm(list=ls())
#------------------------------------
#packages
library(TransPhylo) #version 1.4.5
library(treeio) #version 1.14.1
library(ape) #version 5.4-1 
library(tidyverse) #version 1.3.0
library(lubridate) #version 1.7.9.2
library(coda) #version 0.19-4
library(lattice) #version 0.20-41

#------------------------------------
##Description:
##Reconstruct a transmission tree with the TransPhylo method

#------------------------------------
#------------------------------------
#Transmission scenario
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

#Number of reference tree considered
j <- 1 #tree considered: number out of 30

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes

for (i in 1:nb_scheme){ #sampling schemes
  
  #Phylogenetic tree reconstructed in BEAST2
  tree_beast <- read.beast(paste0("./MCC_",substr(samp,1,1),"/MCC_",samp,"_",j,"_",i,".tree"))
  tree <- as.phylo(tree_beast)
  
  #Sampled sequences
  data <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",samp,"_",j,"_",i,".csv"))
  tstart <- as.Date("01/01/07", "%d/%m/%y") #start of the outbreak
  date <- tstart %m+% months(max(data$tremoved))
  
  #Phylo object for Transphylo
  phylo_tree <- ptreeFromPhylo(tree, dateLastSample = decimal_date(date))
  
  #Parameters for the generation and sampling distributions
  distr_gen <- read.csv(paste0("./samp_",substr(samp,1,1),"/For_transphylo/generation_param_",samp,".csv"))
  distr_gen <- distr_gen %>% filter(sim==j)
  distr_samp <- read.csv(paste0("./samp_",substr(samp,1,1),"/For_transphylo/sampling_param_",samp,".csv"))
  distr_samp <- distr_samp %>% filter(sim==j)
  
  #Implement Transphylo
  set.seed(2)
  
  res<-inferTTree(phylo_tree,mcmcIterations=5e5,
                  thinning= 50, 
                  w.mean=distr_samp$mu,w.std=distr_samp$sd,
                  ws.mean=distr_samp$mu,ws.std=distr_samp$sd,
                  dateT=2020.1)
  
  #Export results
  mcmc <- convertToCoda(res, burnin=0.2) #MCMC chain
  med <- medTTree(res, burnin=0.2) #medoid tree
  tree <- as.data.frame(extractTTree(med)$ttree)
  tree <- tree %>% 
    rename("tinfection"="V1", "tremoved"="V2", "infector_id"="V3") %>%
    mutate("infected_id"=c(extractTTree(med)$nam, rep("unsampled", (nrow(tree)-nrow(data)))))  
  
  #Write output
  #MCMC
  write_csv(as.data.frame(mcmc), paste0("mcmc_",samp,"_",j,"_",i,".csv"))
  
  #Medoid tree
  write_csv(tree, paste0("ttree_",samp,"_",j,"_",i,".csv"))
  
  #Transmission matrix 
  WIW <- as.data.frame(computeMatWIW(res))
  write_csv(WIW, paste0("WIW_",samp,"_",j,"_",i,".csv"))
}

