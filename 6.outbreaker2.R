rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(phylotools)
library(ape)
library(outbreaker2)

#------------------------------------
##Description:
##Reconstruct consensus transmission tree with the outbreaker2 method

#------------------------------------
#------------------------------------
#Transmission scenario
samp <- "B1"

#Number of reference tree considered
j <- 1

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1)

#Generation and sampling distributions
distr_gen <- read_csv(paste0("./samp_",substr(samp,1,1),"/For_outbreaker/generation_time_",samp,"_",j,".csv")) 
distr_gen <- as.vector(distr_gen$x)
distr_samp <- read_csv(paste0("./samp_",substr(samp,1,1),"/For_outbreaker/sampling_time_",samp,"_",j,".csv")) 
distr_samp <- as.vector(distr_samp$x)


for (i in 1:nb_scheme){ #sampling schemes
  data <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",samp,"_",j,"_",i,".csv")) 
  
  #Genetic data
  donnees <- data[,-c(2,3)]
  donnees <- donnees %>% 
    rename(seq.name=infected, seq.text=nucl)
  data_fas <- dat2fasta(donnees, "F7_m.fasta")
  data_gen <- read.dna(file= "F7_m.fasta", format="fasta")
  
  #Create outbreaker_data
  dataO <- outbreaker_data(dna = data_gen, dates = data$tremoved, w_dens = distr_gen, f_dens=distr_samp, ids=data$infected)
  
  #Customising settings and priors
  min_d <- min(data$tremoved)
  config2 <- create_config(n_iter = 1e5, sample_every=50, pb = T, find_import=F, min_date=-min_d, prior_pi=c(1,1))
  
  ##Outbreaker function
  set.seed(2)
  
  res <- outbreaker(data = dataO, config2)
  
  #Consensus tree
  tree <- summary(res, burnin=1e4)$tree
  
  #Parameter values from the MCMC chain
  res <- res[, c(2,5,6)]
  
  #---------------Export results----------------------
  
  #MCMC chain
  write_csv(res, paste0("res_",samp,"_",j,"_",i,".csv")) 
  
  #Consensus tree
  write_csv(tree, paste0("res_tree",samp,"_",j,"_",i,".csv")) 
}


