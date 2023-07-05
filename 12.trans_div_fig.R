rm(list=ls())
#------------------------------------
#package
library(tidyverse)
library(ggplot2)

#------------------------------------
##Description:
##Summarize genetic diversity in simulated data:
##>Proportion of unique sequences
##>Mean transmission divergence
##>Figure with number of SNPs separating transmission pairs

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")


#Get results from trans_divergence code
out <- NULL
for (i in c("A1","B1","B2","S1","S4")){
  samp <- i
  out_t1 <- read_csv(paste0("./Trans_div/trans_div_",samp,".csv"))
  out <- rbind(out,out_t1)
}

#Get results from TransPhylo in order to select only those that converged
trans <- NULL
for (samp in c("A1","B1","B2","S1","S4")){
  out_t1 <- read_csv(paste0("./Accuracy/transphylo_WIW_",samp,"_new_acc.csv"))
  trans <- rbind(trans,out_t1)
}

#Select only those that converged
out <- out %>% mutate(conv=ifelse(tree %in% trans$sim, 1, 0)) %>%
  filter(scenario==1)

out <- out %>% mutate(scen_sim=substr(tree, 1, 2))


#Correct names for transmission scenarios
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "Dead-end host",
                          out$scen_sim == "B1" ~ "Reference",
                          out$scen_sim == "B2" ~ "Badger index",
                          out$scen_sim == "S1" ~ "Single-host",
                          out$scen_sim == "S4" ~ "High mutation rate")

#------------------------------------
#Proportion of unique sequences
tab <- distinct(out %>% select(prop_unique, tree, scen_sim, conv))

tab <- tab %>% group_by(scen_sim, conv) %>%
  summarise(median=round(median(prop_unique),1),
            min=round(min(prop_unique),1),
            max=round(max(prop_unique),1))

#------------------------------------
#Mean transmission divergence
tab2 <- distinct(out %>% select(mean_div, tree, scen_sim, conv))

tab2 <- tab2 %>% group_by(scen_sim, conv) %>%
  summarise(median=round(median(mean_div),2),
            min=round(min(mean_div),2),
            max=round(max(mean_div),2))

#------------------------------------
#Figure in supplementary

out$dist_gen <- as.character(out$dist_gen)

out$scen_sim <- factor(out$scen_sim, levels=c("Reference", "High mutation rate",
                                              "Single-host", "Dead-end host", "Badger index"))

(A <- ggplot(data=out,aes(x=dist_gen, y=nb))+
    geom_boxplot(width=.5)+
    facet_grid(~ scen_sim)+
    theme_bw(base_size=12)+
    labs(x="Genetic distance between transmission pairs", 
    y="Proportion of transmission pairs"))


ggsave("Trans_div.tiff",device="tiff", dpi="print", width=17, height=10, units="cm",  compression="lzw", scale=1)

