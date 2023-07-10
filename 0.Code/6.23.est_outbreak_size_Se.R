rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(broom) #version 0.7.2
library(MASS) #version 7.3-53

#------------------------------------
##Description:
##Alternative transmission scenarios: S1 (single-host) and S4 (high mutation rate)
##Negative Binomial GLMM on the outbreak size indicator
##Comparison with subtree size

#------------------------------------
#------------------------------------
#import tree size results from every method
out <- NULL
for (i in c("S1","S4")){ #S1 (single-host) and S4 (high mutation rate)
  samp <- i
  
  out_o1 <- read_csv(paste0("./Outbreak_size/outbreaker_",samp,"_tree_size.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./Outbreak_size/transphylo_",samp,"_tree_size.csv"))
  out_t1$method <- "transphylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(tree %in% out_t1$tree)
  out <- rbind(out,tot)
}


#scen_sim stands for transmission scenarios
out <- out %>% filter(scenario==1) %>%
  mutate(scen_sim=substr(tree, 1, 2))


#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "S1" ~ "Oc", #meaning only cattle
                          out$scen_sim == "S4" ~ "CTrW_H") #meaning cattle index, wild boars that transmit and high mutation rate

#------------------------------------
#Select transmission scenario 
out <- out %>% filter(scen_sim=="Oc") #either "Oc" or "CTrW_H"

#Implement the GLM on the outbreak size results
model <- glm.nb(nb_sim 
                  ~ offset(log(nb_ref))+ method +0,
                  data=out)

summary(model) #and to get p-value

#Create data.frame with results (S11 Table)
all_meth <- tidy(model, exponentiate = TRUE) %>% 
  mutate(IRR=round(estimate, 2),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3))) 
