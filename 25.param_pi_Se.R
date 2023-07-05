rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(broom)
library(MASS)

#------------------------------------
##Description:
##Negative Binomial GLMM on the outbreak size indicator
##Comparison with reference tree size

#------------------------------------
#------------------------------------
#Estimate number of sequences per tree
nb <- NULL
for (samp in c("A1", "A2", "B1", "B2", "S1","S2","S3","S4","S5","S6")){
  a <- substr(samp, 1, 1)
  for (j in 1:30){
      seq <- read_csv(paste0("./seq_",a,"/seq_",samp,"_",j,"_1.csv"))
      tot <- tibble(nb=nrow(seq), sim=paste(samp,j,sep="_"), scenario="1")
      nb <- rbind(nb,tot)
  }
}


#import parameter results from every method
out <- NULL
for (i in c("A1", "A2", "B1", "B2", "S1","S2","S3","S4","S5","S6")){
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

#Only consider pi parameter and sampling scheme A
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
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "CNTrW",
                          out$scen_sim == "A2" ~ "BNTrW",
                          out$scen_sim == "B1" ~ "CTrW",
                          out$scen_sim == "B2" ~ "BTrW",
                          out$scen_sim == "S1" ~ "Oc",
                          out$scen_sim == "S2" ~ "CNTrW_H",
                          out$scen_sim == "S3" ~ "BNTrW_H",
                          out$scen_sim == "S4" ~ "CTrW_H",
                          out$scen_sim == "S5" ~ "BTrW_H",
                          out$scen_sim == "S6" ~ "Oc_H")
out$scen_sim <- factor(out$scen_sim, levels=c("BNTrW", "BTrW", "CNTrW", "CTrW", "Oc", "BNTrW_H", "BTrW_H", "CNTrW_H", "CTrW_H", "Oc_H"))

#Write output for figure outbreak size
#write_csv(out, "res_outbreak_size_pi_Se.csv")

#------------------------------------
#Select transmission scenario 
out <- out %>% filter(scen_sim=="Oc")

#Implement the GLM on the outbreak size results 
model <- glm.nb(est_size 
                ~ offset(log(real_size))+ method +0,
                data=out)

summary(model) #and to get p-value

#Create data.frame with results
all_meth <- tidy(model, exponentiate = TRUE) %>% 
  mutate(IRR=round(estimate, 2),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3))) 
