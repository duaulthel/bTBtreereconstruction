rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(broom)
library(lme4)

#------------------------------------
##Description:
##Negative Binomial GLMM on the host contribution indicator
##Comparison with reference tree

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")


#import parameter results from every method
samp <- "B1"

out_o1 <- read_csv(paste0("./mcmc/outbreaker_",samp,"_mcmc.csv"))
out_o1$method <- "outbreaker2"

out_t1 <- read_csv(paste0("./mcmc/transphylo_",samp,"_mcmc.csv"))
out_t1$method <- "transphylo"

out <- rbind(out_o1, out_t1)

#keep only trees that converged in TransPhylo
out <- out %>% filter(sim %in% out_t1$sim)


#Only consider pi parameter
pi_param <- out %>% 
  filter(param=="pi") 

class(pi_param$scenario) <- "character"


#import host contribution results from every method
out_o1 <- read_csv(paste0("./Host_contribution/outbreaker_",samp,"_reff_pi.csv"))
out_o1$method <- "outbreaker2"

out_t1 <- read_csv(paste0("./Host_contribution/transphylo_",samp,"_reff_pi.csv"))
out_t1$method <- "transphylo"

out <- rbind(out_o1, out_t1)

#keep only trees that converged in TransPhylo
out <- out %>% filter(tree %in% out_t1$tree)

#Estimate sum
out <- out %>% mutate(R=R*nb, R_sim=R_sim*nb_sim)

out <- out %>% rename(sim=tree) #in order to join out and pi

class(out$scenario) <- "character"

#Join results from pi parameter and host contribution
out <- full_join(out, pi_param, by=c("sim", "scenario", "method"))

#Estimate mean number of secondary cases divided by pi
out <- out %>% mutate(R_sim=round(R_sim/median, 0))

#Correct names for sampling scenarios
out$scenario <- case_when(out$scenario == 1 ~ "A",
                          out$scenario == 2 ~ "T",
                          out$scenario == 4 ~ "SB",
                          out$scenario == 3 ~ "SW",
                          out$scenario == 5 ~ "T+SW",
                          out$scenario == 6 ~ "T+SB")
out$scenario <- factor(out$scenario, levels=c("A", "T", "SB", "T+SB", "SW", "T+SW"))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(sim, 1, 2))

#Correct names for transmission scenarios
out$scen_sim <- "CTrW"

#Random effect variable
out$sim <- as.factor(out$sim)

#Write output for figure outbreak size
write_csv(out, "res_contribution_pi.csv")

#Select host-species and sampling schemes
out <- out %>% filter(sp_infector=="boar") %>%
  filter(scenario!="SW" & scenario!="T+SW") 

#------------------------------------
#Implement the GLMM on the host contribution results
model <- glmer.nb(R_sim 
                  ~ offset(log(R)) + method + method:scenario + (1|sim)+0,
                  control = glmerControl(optimizer="bobyqa"),
                  data=out)

summary(model) #and to get p-value

#Estimate 95% confidence interval
cc <- confint(model,parm="beta_",method="Wald")
ctab <- cbind(est=fixef(model),cc)

#Exponentiate results to get IRR
rtab <- exp(ctab)

#Create data.frame with results
tab <- as.data.frame(round(rtab, 2))

