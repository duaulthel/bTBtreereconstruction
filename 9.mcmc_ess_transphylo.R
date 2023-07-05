rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(coda)

#------------------------------------
##Description:
##Check for convergence and independence of sampling
##& parameter estimation 
##(posterior probability, mutation rate and sampling proportion) 

#------------------------------------
#------------------------------------

#------------------------------------
#Transmission scenario
samp <- "B1"

#Number of sampling schemes
nb_scheme <- ifelse(samp=="B1", 6, 1)

for (j in 1:30){#trees that converged
  
  out_conv <- NULL
  sim <- paste0(samp,"_",j)
  
  for (i in 1:nb_scheme){ #sampling schemes
    
    #results from TransPhylo: mcmc
    res <- read_csv(paste0("./TransPhylo_",samp,"/mcmc_",sim,"_",i,".csv")) 
    
    chain <- mcmc(res) #use coda to transform into mcmc object
    plot(chain) #check trace
    
    ess <- effectiveSize(chain) #estimate ESS
    
    #Create tibble with useful results (median and 95%HPD of parameters)
    est <- tibble("param"=c("pi", "neg", "off.r"), 
                  "median"=c(median(res$pi), median(res$neg), median(res$off.r)),
                  "2.5%"=c(quantile(res$pi, 0.025), quantile(res$neg, 0.025), quantile(res$off.r, 0.025)),
                  "97.5%"=c(quantile(res$pi, 0.975), quantile(res$neg, 0.975), quantile(res$off.r, 0.975)),
                  "ess"=ess[1:3],
                  "sim"=paste(samp, j, sep="_"), 
                  "scenario"=i) #tibble with results for each parameter
    out_conv <- rbind(out_conv, est)
  }
  #Write output for this transmission scenario
  if(j > 1){ #if this is not the first tree, call previous file
    prev <- read_csv(paste0("transphylo_",samp,"_mcmc.csv"))
    out_conv <- rbind(prev, out_conv)
    write_csv(out_conv, paste0("transphylo_",samp,"_mcmc.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(out_conv, paste0("transphylo_",samp,"_mcmc.csv"))
  }
}
