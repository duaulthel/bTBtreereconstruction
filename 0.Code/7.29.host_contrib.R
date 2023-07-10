rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(broom) #version 0.7.2
library(lme4) #version 1.1-29

#------------------------------------
##Description:
##Negative Binomial GLMM on the host contribution indicator
##Comparison with induced subtree 

#------------------------------------
#------------------------------------
#import host contribution results from every method

samp <- "B1" #B1 (reference scenario) 

out_o1 <- read_csv(paste0("./Host_contribution/outbreaker_",samp,"_reff.csv"))
out_o1$method <- "outbreaker2"

out_t1 <- read_csv(paste0("./Host_contribution/transphylo_",samp,"_reff.csv"))
out_t1$method <- "transphylo"

out <- rbind(out_o1, out_t1)

#keep only trees that converged in TransPhylo
out <- out %>% filter(tree %in% out_t1$tree)

out[is.na(out)] <- 0
class(out$scenario) <- "character"

#Correct names for sampling scenarios
out$scenario <- case_when(out$scenario == 1 ~ "Reference",
                          out$scenario == 2 ~ "T",
                          out$scenario == 4 ~ "SB",
                          out$scenario == 3 ~ "SW",
                          out$scenario == 5 ~ "T+SW",
                          out$scenario == 6 ~ "T+SB")
out$scenario <- factor(out$scenario, levels=c("Reference", "T", "SB", "T+SB", "SW", "T+SW"))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(tree, 1, 2))

#Correct names for transmission scenarios
out$scen_sim <- "CTrW" #meaning cattle index and wild boars that transmit

#Random effect variable
out$tree <- as.factor(out$tree)

#Estimate sum
out <- out %>% mutate(R=R*nb, R_sim=R_sim*nb_sim)

#Write output for figure outbreak size
#write_csv(out, "res_contribution.csv")

#Select host-species and sampling schemes
out <- out %>% filter(sp_infector=="boar") %>%
  filter(scenario!="SW" & scenario!="T+SW") 

#------------------------------------
#Implement the GLMM on the host contribution results
model <- glmer.nb(R_sim 
             ~ offset(log(R)) + method + method:scenario + (1|tree)+0,
             control = glmerControl(optimizer="bobyqa"),
             data=out)

summary(model) #and to get p-value

#Estimate 95% confidence interval
cc <- confint(model,parm="beta_",method="Wald")
ctab <- cbind(est=fixef(model),cc)

#Exponentiate results to get IRR
rtab <- exp(ctab)

#Create data.frame with results (S2 Appendix)
tab <- as.data.frame(round(rtab, 2))





