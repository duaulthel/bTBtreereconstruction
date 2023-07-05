rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(ggplot2)

#------------------------------------
##Description:
##Estimate 95%HPD interval for outbreak size
##Comparison with reference tree size

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#Estimate number of sequences per set of cases
nb <- NULL
for (samp in c("B1", "S1", "S4")){
  a <- substr(samp, 1, 1)
  nb_scheme <- ifelse(samp=="B1", 6, 1)
  for (j in 1:30){ #all trees in transmission scenarios
    for (i in 1:nb_scheme){ #all sampling schemes
      seq <- read_csv(paste0("./seq_",a,"/seq_",samp,"_",j,"_",i,".csv"))
      tot <- tibble(nb=nrow(seq), sim=paste(samp,j,sep="_"), scenario=as.character(i))
      nb <- rbind(nb,tot)
    }
  }
}


#import parameter results from every method
out <- NULL
for (samp in c("B1", "S1", "S4")){
  
  out_o1 <- read_csv(paste0("./mcmc/outbreaker_",samp,"_mcmc.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./mcmc/transphylo_",samp,"_mcmc.csv"))
  out_t1$method <- "TransPhylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(sim %in% out_t1$sim)
  out <- rbind(out,tot)
}

#Only consider pi parameter
out <- out %>% 
  filter(param=="pi") 

class(out$scenario) <- "character"

colnames(out) [3] <- "inf"
colnames(out) [4] <- "sup"

#Merge results and proportion of sequences per trees  
out <- left_join(out, nb, by=c("sim", "scenario"))

#Add real outbreak size
out <- out %>% group_by(sim, method) %>% 
  mutate(real_size=nb[scenario==1],
         est_size_sup= round(nb/inf,0),
         est_size_inf= round(nb/sup,0))

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
out$scen_sim <- case_when(out$scen_sim == "B1" ~ "CTrW",
                          out$scen_sim == "S1" ~ "Oc",
                          out$scen_sim == "S4" ~ "CTrW_H")


#Is the real contribution in the 95%HPD?
out <- out %>%
  mutate(correct=case_when(real_size<est_size_inf ~ "Overestimation",
                           real_size>est_size_sup ~ "Underestimation",
                           TRUE ~ "Correct"))

#------------------------------------
#Color code
strans <- c("Underestimation" = "steelblue",
            "Correct" = "chartreuse2",
            "Overestimation" = "indianred3")

#------------------------------------
#Select transmission scenario CTrW only
fig_1 <- out %>% filter(scen_sim=="CTrW")


(A <- ggplot(data=fig_1)+
    geom_errorbar(aes(y=reorder(sim, real_size), xmin=est_size_sup, xmax=est_size_inf, color=correct))+
    geom_point(aes(y=reorder(sim, real_size), x= real_size, color = correct))+
    facet_grid(scenario~method)+
    scale_x_log10()+
    scale_color_manual(name="Estimation", values=strans)+
    theme_bw(base_size=11)+
    labs(x="Number of hosts", y=""))

A + theme(axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          legend.position="bottom")

ggsave("Within_size.tiff",device="tiff", dpi="print", width=17, height=22, units="cm",  compression="lzw", scale=1)

#------------------------------------
#Select alternative transmission scenarios 
fig_2 <- out %>% filter(scen_sim=="CTrW_H") #either CTrW_H or Oc


(A <- ggplot(data=fig_2)+
    geom_errorbar(aes(y=reorder(sim, real_size), xmin=est_size_sup, xmax=est_size_inf, color=correct))+
    geom_point(aes(y=reorder(sim, real_size), x= real_size, color = correct))+
    facet_grid(~method)+
    scale_x_log10()+
    scale_color_manual(name="Estimation", values=strans)+
    theme_bw(base_size=11)+
    labs(x="Number of hosts", y=""))

A + theme(axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          legend.position="bottom")


ggsave("Within_size_high.tiff",device="tiff", dpi="print", width=17, height=9, units="cm",  compression="lzw", scale=1)

