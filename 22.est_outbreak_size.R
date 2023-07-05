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
##Comparison with subtree size

#------------------------------------
#------------------------------------
#import tree size results from every method

samp <- "B1"

out_o1 <- read_csv(paste0("./Outbreak_size/outbreaker_",samp,"_tree_size.csv"))
out_o1$method <- "outbreaker2"

out_t1 <- read_csv(paste0("./Outbreak_size/transphylo_",samp,"_tree_size.csv"))
out_t1$method <- "transphylo"

out <- rbind(out_o1, out_t1)

#keep only trees that converged in TransPhylo
out <- out %>% filter(tree %in% out_t1$tree)

class(out$scenario) <- "character"

#Correct names for sampling scenarios
out$scenario <- case_when(out$scenario == 1 ~ "A",
                          out$scenario == 2 ~ "T",
                          out$scenario == 4 ~ "SB",
                          out$scenario == 3 ~ "SW",
                          out$scenario == 5 ~ "T+SW",
                          out$scenario == 6 ~ "T+SB")
out$scenario <- factor(out$scenario, levels=c("A", "T", "SB", "T+SB", "SW", "T+SW"))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(tree, 1, 2))

#Correct names for transmission scenarios
out$scen_sim <- "CTrW"

#Write output for figure outbreak size
write_csv(out, "res_outbreak_size.csv")

#Random effect variable
out$tree <- as.factor(out$tree)

#------------------------------------
#Implement the GLMM on the outbreak size results
model <- glmer.nb(nb_sim 
             ~ offset(log(nb_ref))+method + method:scenario + (1|tree)+0,
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
