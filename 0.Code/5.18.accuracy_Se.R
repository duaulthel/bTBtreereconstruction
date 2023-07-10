rm(list=ls())
#------------------------------------
#packages
library(tidyverse) #version 1.3.0
library(broom) #version 0.7.2

#------------------------------------
##Description:
##Alternative transmission scenarios: S1 (single-host) and S4 (high mutation rate)
##Binomial GLM on the accuracy indicator

#------------------------------------
#------------------------------------
#import accuracy results from every method
out <- NULL
for (i in c("S1","S4")){ #S1 (single-host) and S4 (high mutation rate)
  samp <- i
  out_s1 <- read_csv(paste0("./Accuracy/seqTrack_",samp,"_acc.csv"))
  out_s1$method <- "seqTrack"
  
  out_o1 <- read_csv(paste0("./Accuracy/outbreaker_",samp,"_acc.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t2 <- read_csv(paste0("./Accuracy/transphylo_WIW_",samp,"_new_acc.csv"))
  out_t2$method <- "TransPhylo"
  
  tot <- rbind(out_s1, out_o1, out_t2)
  out <- rbind(out,tot)
}

#keep only trees that converged in TransPhylo
class(out$scenario) <- "character"
indice <- out %>% filter(method=="TransPhylo")
out <- out %>% filter(sim %in% indice$sim) %>%
  filter(scenario==1) #only keep reference scheme

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(sim, 1, 2))

#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "S1" ~ "Oc", #meaning only cattle
                          out$scen_sim == "S4" ~ "CTrW_H") #meaning cattle index, wild boars that transmit and high mutation rate

#------------------------------------
#Table with results from all methods (S8 Table)
tab <- out %>% 
  mutate(prop=nb_acc/(nb_acc+nb_non_acc)*100)

tab <- tab %>% group_by(scen_sim, method) %>%
  summarise(median=round(median(prop),1),
            min=round(min(prop),1),
            max=round(max(prop),1))

tab <- tab %>% mutate(res=paste0(median, " (", min, "-", max, ")")) %>%
  select(scen_sim, method, res)

tab <- tab %>% pivot_wider(names_from=scen_sim, values_from=res)

write_csv(tab, "res.csv")
#------------------------------------
#Select transmission scenario 
out <- out %>% filter(scen_sim=="Oc") #either "Oc" or "CTrW_H"

#Implement the GLM on the accuracy results
model <- glm(cbind(nb_acc, nb_non_acc) 
             ~ method,
             data=out,
             family=binomial)

summary(model)

#Calculate 95% confidence interval and exponentiate to get Odds Ratio (S9 Table)
all_meth <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(OR=round(estimate, 2),
         IC=paste(round(conf.low,2), round(conf.high,2), sep="-"),
         p.value=ifelse(p.value<0.001, "<0.001", round(p.value,3))) %>%
  select(term, OR, IC, p.value)



