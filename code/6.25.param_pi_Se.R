rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(broom)
library(MASS)

#------------------------------------
##Description: 
##Alternative transmission scenarios
##Negative Binomial GLMM on the outbreak size indicator
##Comparison with reference tree size

#------------------------------------
#------------------------------------
#Estimate number of sequences per tree
nb <- NULL
for (samp in c("S1","S4")){
  a <- substr(samp, 1, 1)
  for (j in 1:30){
      seq <- read_csv(paste0("./seq_",a,"/seq_",samp,"_",j,"_1.csv"))
      tot <- tibble(nb=nrow(seq), sim=paste(samp,j,sep="_"), scenario="1")
      nb <- rbind(nb,tot)
  }
}


#import parameter results from every method
out <- NULL
for (i in c("S1","S4")){
  samp <- i
  
  out_o1 <- read_csv(paste0("./mcmc/outbreaker_",samp,"_mcmc.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./mcmc/transphylo_",samp,"_mcmc.csv"))
  out_t1$method <- "transphylo"
  # 
  tot <- rbind(out_o1, out_t1)
  tot <- tot %>% filter(sim %in% out_t1$sim)
  out <- rbind(out,tot)
}

out$scenario <- as.character(out$scenario)

#Only consider pi parameter and reference sampling scheme 
out <- out %>% 
  filter(scenario==1) %>%
  filter(param=="pi") 

#Merge results and proportion of sequences per trees 
out <- left_join(out, nb, by=c("sim", "scenario"))

#Add real outbreak size
out <- out %>% group_by(sim, method) %>% 
  mutate(real_size=nb[scenario=="1"],
         est_size= round(nb/median,0))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(sim, 1, 2))


#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "S1" ~ "Oc",
                          out$scen_sim == "S4" ~ "CTrW_H")
                          
#Write output for figure outbreak size
#write_csv(out, "res_outbreak_size_pi_Se.csv")

#------------------------------------
#Select transmission scenario 
out <- out %>% filter(scen_sim=="Oc") #Oc or CTrW_H

#Implement the GLM on the outbreak size results 
model <- glm.nb(est_size 
                ~ offset(log(real_size))+ method +0,
                data=out)

summary(model) #and to get p-value

#Create data.frame with results
all_meth <- tidy(model, exponentiate = TRUE) %>% 
  mutate(IRR=round(estimate, 2),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3))) 
