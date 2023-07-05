rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(ggplot2)

#------------------------------------
##Description:
##Estimate number of trees with super-spreaders
##Figure with maximum number of transmissions due to a super-spreader

#------------------------------------
#------------------------------------
#import results from every method
out <- NULL
for (i in c("A1","B1","B2","S1","S4")){
  samp <- i
  out_s1 <- read_csv(paste0("./Super-spreader/seqTrack_",samp,"_ssp.csv"))
  out_s1$method <- "seqTrack"
  
  out_o1 <- read_csv(paste0("./Super-spreader/outbreaker_",samp,"_ssp.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./Super-spreader/transphylo_",samp,"_ssp.csv"))
  out_t1$method <- "TransPhylo"

  tot <- rbind(out_s1, out_o1, out_t1)
  class(out$scenario) <- "character"
  out <- rbind(out,tot)
}

#keep only trees that converged in TransPhylo and 
#that had superspreaders
indice <- out %>% filter(method=="TransPhylo")
out <- out %>% filter(sim %in% indice$sim) %>% 
  filter(!is.na(infector)) %>% #remove trees w/o superspreaders
  mutate(prop=n/nub_trans) #proportion of transmission events concerned

#Correct names for sampling scenarios
out$scenario <- case_when(out$scenario == 1 ~ "Reference",
                          out$scenario == 2 ~ "T",
                          out$scenario == 4 ~ "SB",
                          out$scenario == 3 ~ "SW",
                          out$scenario == 5 ~ "T+SW",
                          out$scenario == 6 ~ "T+SB")
out$scenario <- factor(out$scenario, levels=c("Reference", "T", "SB", "T+SB", "SW", "T+SW"))

#scen_sim stands for transmission scenarios
out <- out %>% mutate(scen_sim=substr(sim, 1, 2))

#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "Dead-end host",
                          out$scen_sim == "B1" ~ "Reference",
                          out$scen_sim == "B2" ~ "Badger index",
                          out$scen_sim == "S1" ~ "Single-host",
                          out$scen_sim == "S4" ~ "High mutation rate")

#--------------------------------------------------------
#Create table with proportion and number of trees with ssp per transmission
#and sampling scenario
tab <- distinct(out %>% select(scen_sim, sim, scenario, method)
                %>% mutate(n = 1))

#Table with the wanted numbers
tab <- tab %>% group_by(method, scenario,scen_sim) %>%
  summarise(nb=sum(n)) 

#---------------------------------------------------
#Only keep the maximum number of secondary cases per super-spreader
tab2 <- out %>% group_by(sim, scenario, scen_sim, method) %>%
  slice_max(n, n= 1) %>% mutate(nb=1)

#Number of trees with a cattle/badger/boar as the most prolific super-spreader
tab3 <- tab2 %>% group_by(scenario, scen_sim, method, sp_infector) %>%
  summarise(prop_sp=sum(nb))

#Median and range of the maximum number of cases due to a single super-spreader
tab2 <- tab2 %>% group_by(scenario, scen_sim, method) %>%
  summarise(median=median(n),
            min=min(n),
            max=max(n))
#---------------------------------------------------

