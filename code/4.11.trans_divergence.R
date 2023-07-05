rm(list=ls())
#------------------------------------
#package
library(ape)
library(phylotools)
library(lubridate)
library(tidyverse)
library(adegenet)
library(reshape2)

#------------------------------------
##Description:
##Estimate genetic diversity in simulated data:
##>Proportion of unique sequences
##>Mean transmission divergence

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1"

for (j in 1:30){#all trees in the transmission scenario
  sim <- paste0(samp,"_",j) #name of reference tree
  
  #reference tree
  Ttree <- read_csv(paste0("./Ttrees_",substr(samp,1,1),"/Ttree_",sim,".csv"))
  
  #remove index case (infector == NA) in order to consider transmission events
  Ttree <- Ttree %>% filter(!is.na(infector))
  
  #Create column with transmission events
  Ttree <- Ttree %>% 
    mutate(trans=paste(infector, infected, sep="->"))
  
  #genetic sequences
  data <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",sim,"_1.csv")) 
  
  #Estimate number of unique sequences
  unique <- data %>% group_by(nucl) %>%
    summarise(nb=n_distinct(infected)) %>% filter(nb==1)
  
  #Create DNA bin object from sequences
  donnees <- data %>% select(infected, nucl)
  donnees <- donnees %>% 
    rename(seq.name=infected, seq.text=nucl) #create fasta file
  data_fas <- dat2fasta(donnees, "F7_m.fasta") #export fasta file
  data_gen <- read.dna(file= "F7_m.fasta", format="fasta") #import fasta file
  
  #Estimate pairwise genetic distances in number of SNPs
  dm_model <- dist.dna(data_gen, model="N", as.matrix=T)
  
  #Attibute distance to piar of infected hosts
  dist_tree <- setNames(melt(as.matrix(dm_model)), c('infector', 'infected', 'dist_gen'))
  
  #Only consider direct transmission pairs
  dist_tree <- dist_tree %>% 
    mutate(trans=paste(infector, infected, sep="->")) %>%
    filter(trans %in% Ttree$trans)
  
  #Create output
  tab_dist <- dist_tree %>% group_by(dist_gen) %>%
    summarise(nb=n_distinct(infected)/nrow(dist_tree)*100,
              prop_unique=nrow(unique)/nrow(dist_tree)*100,
              tree=sim,
              scenario=i)
  
  #Add mean transmission divergence
  tab_dist <- cbind(tab, tibble("mean_div"=mean(dist_tree$dist_gen)))

 
  #-----------------------OUTPUT-----------------------------------
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("trans_div_",samp,".csv"))
    tab_dist <- rbind(prev, tab_dist)
    write_csv(tab_dist, paste0("trans_div_",samp,".csv"))
  }else{ #if this is the first tree, write new file
    write_csv(tab_dist, paste0("trans_div_",samp,".csv"))
  }
}

