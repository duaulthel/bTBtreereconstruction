rm(list=ls())
#------------------------------------
#packages
library(ape) #version 5.4-1
library(phylotools) #version 0.2.2
library(lubridate) #version 1.7.9.2
library(tidyverse) #version 1.3.0
library(adegenet) #version 2.1.3

#------------------------------------
##Description:
##Reconstruct a transmission tree with the seqTrack method

#------------------------------------
#------------------------------------
#start of the outbreak
tstart <- as.Date("01/01/07", "%d/%m/%y")

#Transmission scenario
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes

for (j in 1:30){ #all trees in the transmission scenario
  sim <- paste0(samp,"_",j) 
  
  for (i in 1:nb_scheme){ #sampling schemes
    
    #genetic sequences
    data <- read_csv(paste0("./seq_",substr(samp,1,1),"/seq_",samp,"_",sim,"_",i,".csv")) 
    data$tremoved <- tstart %m+% months(data$tremoved)
    
    #Pairwise genetic distances
    donnees <- data %>% select(infected, nucl)
    donnees <- donnees %>% 
      rename(seq.name=infected, seq.text=nucl)
    data_fas <- dat2fasta(donnees, "F7_m.fasta")
    data_gen <- read.dna(file= "F7_m.fasta", format="fasta")
    dm_model <- dist.dna(data_gen, model="F84", as.matrix=T) #substitution model= F84
    
    #Implement seqtrack
    res <- seqTrack(dm_model, x.names=data$infected, x.dates=data$tremoved) 
    
    #Create transmission tree output with information on hosts
    res <- res %>% rownames_to_column(var="infected")
    get_res <- left_join(res, data[,c(1,3)], by="infected")
    
    get_res <- get_res %>% mutate("sp_infector"=sp[ances],
                                  "infector"=infected[ances])
    
    #------------------------------------
    #Write output
    write_csv(get_res, paste0("seqtrack_tree",sim,"_",i,".csv"))
  }
  
}



