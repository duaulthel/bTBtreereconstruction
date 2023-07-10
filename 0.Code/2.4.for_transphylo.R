rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(epitrix) #version 0.2.2
library(lubridate) #version 1.7.9.2
library(fitdistrplus) #version 1.1-6

#------------------------------------
##Description:
##Estimates shape, mean and sd for gamma distributions
##of generation times and sampling intervals needed to implement TransPhylo 

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "A1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

tstart <- as.Date("01/01/07", "%d/%m/%y") #start of the outbreak

distr_gen <- NULL
distr_samp <- NULL


for (j in 1:30){ #all trees in transmission scenario
  #------------------------------------
  #Generation time distribution
  
  #reference tree
  Ttree <- read_csv(paste0("./Ttrees_",substr(samp,1,1),"/Ttree_",samp,"_",j,".csv")) 
  
  #For Transphylo: change months into years
  Ttree$tinfection <- tstart %m+% months(Ttree$tinfection)
  Ttree$tinfection <-  decimal_date(Ttree$tinfection)
  
  #Get generation times in the reference tree
  out <- NULL
  for(l in Ttree$infected[-1]){
    sub_tree <- Ttree %>% filter(infected == l) 
    gen_inf <- sub_tree$infector  
    gen_time <- sub_tree$tinfection-Ttree[Ttree$infected == gen_inf,]$tinfection
    out <- c(out, gen_time)
  }
  
  #Fit a gamma distribution to the generation times
  si_fit <- fit_disc_gamma(out)
  
  #Output with parameters
  distr_gen <- rbind(distr_gen, tibble("mu"=si_fit$mu, 
                                       "sd"=si_fit$sd, 
                                       "shape"=si_fit$distribution$parameters$shape, 
                                       "sim"=j))
  
  #Sampling time distribution
  
  #For Transphylo: change months into years
  Ttree$tremoved <- tstart %m+% months(Ttree$tremoved)
  Ttree$tremoved <-  decimal_date(Ttree$tremoved)
  
  #Get sampling intervals in the reference trees
  out <- NULL
  for(k in Ttree$infected){
    sub_tree <- Ttree %>% filter(infected == k) 
    samp_time <- sub_tree$tremoved-sub_tree$tinfection
    out <- c(out, samp_time)
  }
  
  #Fit a gamma distribution
  si_fit <- fit_disc_gamma(out)
  
  #Output with parameters
  distr_samp <- rbind(distr_samp, tibble("mu"=si_fit$mu, 
                                         "sd"=si_fit$sd, 
                                         "shape"=si_fit$distribution$parameters$shape, 
                                         "sim"=j))

}

#

#------------------------------------
#Write output

#Generation times
write.csv(distr_gen, paste0("generation_param_",samp,".csv"))

#Sampling intervals
write.csv(distr_samp, paste0("sampling_param_",samp,".csv")) 

