rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(broom) #version 0.7.2
library(MASS) #version 7.3-53

#------------------------------------
##Description:
##Alternative transmission scenarios: A1 (dead-end) and S4 (high mutation rate)
##Negative Binomial GLMM on the host contribution indicator
##Comparison with reference tree

#------------------------------------
#------------------------------------
#import parameter results from every method
out <- NULL
for (i in c("A1", "S4")){ #A1 (dead-end) and S4 (high mutation rate)
  samp <- i
  
  out_o1 <- read_csv(paste0("./mcmc/outbreaker_",samp,"_mcmc.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./mcmc/transphylo_",samp,"_mcmc.csv"))
  out_t1$method <- "transphylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(sim %in% out_t1$sim)
  out <- rbind(out,tot)
}

#Only consider pi parameter
pi_param <- out %>% 
  filter(param=="pi") %>%
  filter(scenario==1)


#import host contribution results from every method
out <- NULL
for (i in c("A1", "S4")){ #A1 (dead-end) and S4 (high mutation rate)
  samp <- i
  
  out_o1 <- read_csv(paste0("./Host_contribution/outbreaker_",samp,"_reff_pi.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./Host_contribution/transphylo_",samp,"_reff_pi.csv"))
  out_t1$method <- "transphylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(tree %in% out_t1$tree)
  out <- rbind(out,tot)
}

out <- out %>% filter(scenario==1) %>% #only keep reference scheme
  rename(sim=tree) #in order to join out and pi

#Estimate sum
out <- out %>% mutate(R=R*nb, R_sim=R_sim*nb_sim)

#Join results from pi parameter and host contribution
out <- full_join(out, pi_param, by=c("sim", "method"))

#Estimate mean number of secondary cases divided by pi
out <- out %>% mutate(R_sim=round(R_sim/median, 0))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(sim, 1, 2))


#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "CNTrW", #meaning cattle index and wild boars that did NOT transmit
                          out$scen_sim == "S4" ~ "CTrW_H") #meaning cattle index, wild boars that transmit and high rate

#Write output for supplementary
write_csv(out, "res_contribution_pi_Se.csv")

#------------------------------------
#Select host-species and transmission scenario 
out <- out %>% filter(sp_infector=="cattle") %>%
  filter(scen_sim=="CNTrW") #either "CNTrW" or "CTrW_H"

#Implement the LM on the host contribution results
model <- glm.nb(R_sim 
              ~ offset(log(R))+method+0,
              data=out)

summary(model) #and to get p-value

#Create data.frame with results (S12 Table)
all_meth <- tidy(model, exponentiate = TRUE) %>%
  mutate(OR=round(estimate, 2),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3)))

