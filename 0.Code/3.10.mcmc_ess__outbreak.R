rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(coda) #version 0.19-4

#------------------------------------
##Description:
##Check for convergence and independence of sampling
##& parameter estimation 
##(posterior probability, mutation rate and sampling proportion) 

#------------------------------------
#------------------------------------

#Transmission scenario considered
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1) #do not change B1, only scenario with 6 schemes

for (j in 1:30){ #all trees in transmission scenario
  
  out_conv <- NULL
  sim <- paste0(samp,"_",j)
  
  for (i in 1:nb_scheme){ #sampling schemes
    
    #results from outbreaker2: mcmc
    res <- read_csv(paste0("./outbreaker_",samp,"/res_",sim,"_",i,".csv")) 
    res <- res[-c(1:round(0.1*nrow(res))),] #remove the 10% burn-in
    
    chain <- mcmc(res) #use coda to transform into mcmc object
    plot(chain) #check trace
    
    ess <- effectiveSize(chain) #estimate ESS
    
    #Create tibble with useful results (median and 95%HPD of parameters)
    est <- tibble("param"=c("post", "mu", "pi"), 
                  "median"=c(median(res$post), median(res$mu), median(res$pi)),
                  "2.5%"=c(quantile(res$post, 0.025), quantile(res$mu, 0.025), quantile(res$pi, 0.025)),
                  "97.5%"=c(quantile(res$post, 0.975), quantile(res$mu, 0.975), quantile(res$pi, 0.975)),
                  "ess"=ess,
                  "sim"=paste(samp, j, sep="_"), 
                  "scenario"=i) #tibble with results for each parameter
    out_conv <- rbind(out_conv, est)
  }
  #Write output for this transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("outbreaker_",samp,"_mcmc.csv"))
    out_conv <- rbind(prev, out_conv)
    write_csv(out_conv, paste0("outbreaker_",samp,"_mcmc.csv"))
  }
  else{ #if this is the first tree, write new file
    write_csv(out_conv, paste0("outbreaker_",samp,"_mcmc.csv"))
  }
}


