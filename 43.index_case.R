rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(broom)
library(lme4)

#------------------------------------
##Description:
##Binomial GLMM on the index case indicator

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#import index case results from every method
out <- NULL
for (i in c("A1","B1")){
  samp <- i
  out_s1 <- read_csv(paste0("./Index_case/seqTrack_",samp,"_index.csv"))
  out_s1$method <- "seqTrack"

  out_o1 <- read_csv(paste0("./Index_case/outbreaker_",samp,"_index.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./Index_case/transphylo_",samp,"_index.csv"))
  out_t1$method <- "TransPhylo"
  
  tot <- rbind(out_s1, out_o1, out_t1)
  out <- rbind(out,tot)
}

#keep only trees that converged in TransPhylo
indice <- out %>% filter(method=="TransPhylo")
out <- out %>% filter(sim %in% indice$sim) 

class(out$scenario) <- "character"

#Transform index outcome into T/F
out$index <- ifelse(out$index==1, T, F)

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
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "CNTrW",
                          out$scen_sim == "B1" ~ "CTrW")

#Random effect variable
out$sim <- as.factor(out$sim)

#------------------------------------
#Table with results from all methods
tab <- out %>% 
  mutate(index=ifelse(index==T, 1, 0),
         n=1)

tab <- tab %>% group_by(scen_sim, scenario, method) %>%
  summarise(prop=round(sum(index)/sum(n)*100, 1))

tab <- tab %>% pivot_wider(names_from=scen_sim, values_from=prop)
#write_csv(tab, "res.csv")

#------------------------------------
#------------------------------------
#Select transmission scenario CTrW only
test <- out %>% filter(scen_sim=="CTrW") 

#Implement the GLM on the host contribution results
model <- glmer(index 
             ~ method + method:scenario +(1|sim),
             data=test,
             family=binomial)

# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model is nearly unidentifiable: large eigenvalue ratio
#                - Rescale variables?


summary(model) #and to get p-value

#Estimate 95% confidence interval
cc <- confint(model,parm="beta_",method="Wald")
ctab <- cbind(est=fixef(model),cc)

#Exponentiate results to get Odds Ratio
rtab <- exp(ctab)

#Create data.frame with results
tab <- as.data.frame(round(rtab, 2))

#------------------------------------
#Select transmission scenario CTrW only
test <- out %>% filter(scen_sim=="CNTrW") 

#Implement the GLM on the host contribution results
model <- glm(index 
               ~ method,
               data=test,
               family=binomial)

summary(model) #and to get p-value

#Calculate 95% confidence interval and exponentiate to get Odds Ratio
all_meth <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(OR=round(estimate, 2),
         IC=paste(round(conf.low,2), round(conf.high,2), sep="-"),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3)))
