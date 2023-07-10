rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(broom) #version 0.7.2
library(MASS) #version 7.3-53

#------------------------------------
##Description:
##Alternative transmission scenarios: A1 (dead-end) and S4 (high mutation rate)
##Negative Binomial GLM on the host contribution indicator
##Comparison with induced subtree 

#------------------------------------
#------------------------------------
#import host contribution results from every method
out <- NULL
for (i in c("A1", "S4")){ #A1 (dead-end) and S4 (high mutation rate)
  samp <- i
  
  out_o1 <- read_csv(paste0("./Host_contribution/outbreaker_",samp,"_reff.csv"))
  out_o1$method <- "outbreaker2"

  out_t1 <- read_csv(paste0("./Host_contribution/transphylo_",samp,"_reff.csv"))
  out_t1$method <- "transphylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(tree %in% out_t1$tree)
  out <- rbind(out,tot)
}


#scen_sim stands for transmission scenarios
out <- out %>% filter(scenario==1) %>% #only keep reference scheme
  mutate(scen_sim=substr(tree, 1, 2))


#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "CNTrW", #meaning cattle index and wild boars that did NOT transmit
                          out$scen_sim == "S4" ~ "CTrW_H") #meaning cattle index, wild boars that transmit and high rate

#Estimate sum
out <- out %>% mutate(R=R*nb, R_sim=R_sim*nb_sim)

#Write output for supplementary
#write_csv(out, "res_contribution_Se.csv")

#------------------------------------
#Select host-species and transmission scenario 
out <- out %>% filter(sp_infector=="badger") %>%
  filter(scen_sim=="CNTrW") #either "CNTrW" or "CTrW_H"

#Implement the LM on the host contribution results
model <- glm.nb(R_sim 
              ~ offset(log(R))+method+0,
              data=out)

summary(model) #and to get p-value

#Create data.frame with results (S12 Table)
all_meth <- tidy(model, exponentiate = TRUE) %>%
  mutate(IRR=round(estimate, 2),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3)))

 
