rm(list=ls())
#------------------------------------
#package
library(tidyverse) #version 1.3.0

#------------------------------------
##Description:
##1. Get parameters needed to simulate sequences.
##2. Simulate sequences along the transmission tree 
##NB: no within-host diversity considered here

#------------------------------------
#------------------------------------

#1. Get parameters needed to simulate sequences.

#------------------------------------
#Transmission scenario considered
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

j <- 1 #tree considered: number out of 30

#------------------------------------
#import transmission tree
Ttree <- read_csv(paste0("./Ttrees_",substr(samp,1,1),"/Ttree_",samp,"_",j,".csv")) %>% 
  mutate(infector=ifelse(is.na(infector), infected, infector))%>% 
  distinct()

#give hypothetical sampling date to infected hosts not yet removed
Ttree[Ttree$tremoved==-1,]$tremoved <- 157 

#Transmission pathways: needed to obtain order of infection and later simulate strains
#function that stores transmission pathways
find_path <- function(res, Ttree){
  out <- res
  while(Ttree[Ttree$infected== res,]$infector != res){
    res <- Ttree[Ttree$infected== res,]$infector
    out <- rbind(res, out)
  }
  return(out)
}

#apply find_path at whole Ttree
Tpath <- NULL
for (i in Ttree$infected){
  Tpath <- rbind(Tpath, find_path(i, Ttree))
}

#weed out duplicates in Tpath and thus obtain order of transmission
Tpath <- unique(Tpath)

#Hosts listed according to time of infection
init_inf <- Tpath[1] #initially infected host
init_suscep <- Tpath[-1] #initially susceptible hosts


#------------------------------------
#Initializing sequences in the region

data <- read_csv2("SNP_167souches_AF2122.csv") #sequences sampled in PA Landes
data <- data[, -c(1,57, 92, 127)] #remove SNPs in repetitive areas of the genome
n_snp <- length(data) # number of SNP


#create data.frame with the sequence infecting the index case
data <- cbind("host_id"=init_inf, 
              "sp"=Ttree[Ttree$infected==init_inf,]$sp_infected, 
              "seq"=paste0(init_inf, "_", 1),
              "time"= 0, 
              "time_flag"="tinf", 
              "modeinf"="init",
              data[113,]) #sample(1:nrow(data), 1) = 113 sampled sequence used for all simulations

#------------------------------------
#Create substitution model

mu.rate <- 0.024*n_snp/12 # mutation rate per genome per month 
#low mutation rate: 0.0024
#high mutation rate: 0.024

##Substitution probabilities
#parameters obtained from SNP_167souches_AF2122.csv 
pi_A= 0.233 
pi_C= 0.274
pi_G= 0.275
pi_T= 0.218
kappa= 5.9 

sum_line1= pi_C+kappa*pi_G+pi_T
sum_line2= pi_A+pi_G+kappa*pi_T
sum_line3= kappa*pi_A+pi_C+pi_T
sum_line4= pi_A+kappa*pi_C+pi_G

p_subst = data.frame("A" = c(0, pi_C/sum_line1, kappa*pi_G/sum_line1, pi_T/sum_line1),
                     "C" = c(pi_A/sum_line2, 0, pi_G/sum_line2, kappa*pi_T/sum_line2),
                     "G" = c(kappa*pi_A/sum_line3, pi_C/sum_line3, 0, pi_T/sum_line3),
                     "T" = c(pi_A/sum_line4, kappa*pi_C/sum_line4, pi_G/sum_line4, 0),
                     row.names = c("A", "C", "G", "T"))
p_subst <- as.matrix(p_subst)

#Possible variation between species 
# (not used in article, hence the =1)
#(coefficient that will be multiplied by mu.rate)
m_c <- 1 #cattle
m_b <- 1 #badger
m_w <- 1 #wild boar


#------------------------------------
#------------------------------------
#2. Simulate sequences along the transmission tree 
# (reminder: no within-host diversity)

#------------------------------------
#First function: mutat = mutates a single SNP
mutat <- function(res, tend) {
  coeff <- case_when(res$sp=="cattle" ~ m_c,  
                     res$sp=="badger" ~ m_b,
                     res$sp=="boar"~ m_w) 
  #coefficient to vary mu.rate according to host species (here = 1)
  
  t_mu <- rexp(1,1/3*mu.rate*coeff) #sample time to mutation
  
  if(res$time + t_mu < tend){ #if time sampled before removal
    res$time <- res$time + t_mu #update time
    idx <- sample(n_snp, size = 1, replace = FALSE) #sample nucleotide position of mutation
    res[,idx+6] <- sample(colnames(p_subst), size=1, prob=p_subst[rownames(p_subst)==as.character(res[,idx+6]),]) 
    #substitute nucleotide in new sequence (+6 since time+host+flag+seq+sp+mode was added)
    
    return(res) #function returns new sequence
    
  }
  else{
    
    return(NULL) #no mutation possible, mutation after host removal
  }
}

#Second function gen.seq= generate sequences during infection in a host
gen.seq <- function(host, Ttree, data) {
  Ttreehost <- Ttree %>% filter(infected==host) #select infected host considered 
  tend <- Ttree[Ttree$infected==host,]$tremoved #find tend to use mutat
  
  poss_inf <- data %>% 
    filter(host_id == Ttreehost$infector) %>% 
    filter(time <= Ttreehost$tinfection) # find possible transmitted strains
  poss_inf <- tail(poss_inf, 1) #take only the most recent strain (no within-host diversity)
  seq_host <- tibble(host_id=host,
                     sp=Ttreehost$sp_infected,
                     seq=paste0(host, "_", 0),
                     time=Ttreehost$tinfection , 
                     time_flag="tinf",
                     modeinf=Ttreehost$modeinfection,
                     poss_inf[,-c(1:6)])
  while (nrow(seq_host %>% filter(time_flag != "end"))!=0){
    seq_poss <- tail(seq_host,1)
    #select strains that can mutate (only the most recent here)
    seq_add <- mutat(seq_poss, tend) #mutate sequence within the host
    
    if(is.null(seq_add)) {seq_poss$time_flag <- "end" #no mutation possible
    } 
    
    seq_host <- rbind(seq_host, seq_add)
    seq_host$seq <- paste(host, 1:nrow(seq_host), sep="_")
    return(seq_host)
  }
  
}

#------------------------------------
#simulate sequences in initially infected hosts of the transmission tree
#separate from other hosts in order not to duplicate data

seq <- NULL
for (i in init_inf){ 
  seq <- rbind(seq, gen.seq(i, Ttree, data))
}

#simulate sequences in initially susceptible hosts of the transmission tree
for (k in init_suscep){
  seq <- rbind(seq, gen.seq(k, Ttree, seq))
}


#------------------------------------
#Sampled sequences from simulated sequences
#One sequence sampled per animal, the most recent (no within-host diversity here)
out <- NULL
for (i in unique(Tpath)){
  sub_seq <- tail(seq[seq$host_id==i,],1)
  out <- rbind(out, sub_seq)
}


out <- out %>% unite(col="nucl", 7:177, sep="") #unite nucleotides

#Replace time of most recent sequence by sampling time (= removal)
out <- right_join(Ttree[,c(6,13)], out[, -c(3,4,5,6)], by=c("infected"="host_id"))



#------------------------------------

write_csv(out, paste0("seq_",samp,"_",j,"_1.csv"))

