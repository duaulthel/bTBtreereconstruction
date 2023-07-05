rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(broom)
library(MASS)
library(lme4)
#------------------------------------
##Description:
##Negative Binomial GLMM on the outbreak size indicator
##Comparison with reference tree size

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#Estimate number of sequences per set of cases

samp <- "B1"
a <- substr(samp, 1, 1)
for (j in 1:30){ #all trees in transmission scenarios
  for (i in 1:6){ #all sampling schemes
    seq <- read_csv(paste0("./seq_",a,"/seq_",samp,"_",j,"_",i,".csv"))
    tot <- tibble(nb=nrow(seq), sim=paste(samp,j,sep="_"), scenario=as.character(i))
    nb <- rbind(nb,tot)
  }
}



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
out <- out %>% 
  filter(param=="pi") 

class(out$scenario) <- "character"


#Merge results and proportion of sequences per trees  
out <- left_join(out, nb, by=c("sim", "scenario"))

#Add real outbreak size
out <- out %>% group_by(sim, method) %>% 
  mutate(real_size=nb[scenario==1],
         est_size= round(nb/median,0))

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

#Write output for figure outbreak size
write_csv(out, "res_outbreak_size_pi.csv")

#Random effect variable
out$sim <- as.factor(out$sim)

#------------------------------------
#Implement the GLMM on the outbreak size results
model <- glmer.nb(est_size 
                ~ offset(log(real_size))+ method+ method:scenario + (1|sim)+0,
                control = glmerControl(optimizer="bobyqa"),
                data=out)

#Continue with previous model fit
ss <- getME(model,c("theta","fixef"))
model2 <- update(model,
                 start=ss,
                 control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e4)))

summary(model2) #and to get p-value

#Estimate 95% confidence interval
cc <- confint(model2,parm="beta_",method="Wald")
ctab <- cbind(est=fixef(model2),cc)

#Exponentiate results to get IRR
rtab <- exp(ctab)

#Create data.frame with results
tab <- as.data.frame(round(rtab, 2))
