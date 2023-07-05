rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(ggplot2)
library(tidytext)

#------------------------------------
##Description:
##Estimate HPD interval
##Comparison with reference tree

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")


#import parameter results from every method
out <- NULL
for (i in c("B1", "A1", "S4")){
  samp <- i
  
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
pi_param <- out %>% 
  filter(param=="pi") 

class(pi_param$scenario) <- "character"

colnames(pi_param) [3] <- "inf"
colnames(pi_param) [4] <- "sup"

#import host contribution results from every method
out <- NULL
for (i in c("B1", "A1", "S4")){
  samp <- i
  
  out_o1 <- read_csv(paste0("./Host_contribution/outbreaker_",samp,"_reff_pi.csv"))
  out_o1$method <- "outbreaker2"
  
  out_t1 <- read_csv(paste0("./Host_contribution/transphylo_",samp,"_reff_pi.csv"))
  out_t1$method <- "TransPhylo"
  
  tot <- rbind(out_o1, out_t1)
  
  #keep only trees that converged in TransPhylo
  tot <- tot %>% filter(tree %in% out_t1$tree)
  out <- rbind(out,tot)
}


#Estimate sum
out <- out %>% 
  mutate(R=R*nb, R_sim=R_sim*nb_sim)

out <- out %>% rename(sim=tree) #in order to join out and pi

class(out$scenario) <- "character"

#Join results from pi parameter and host contribution
out <- full_join(out, pi_param, by=c("sim", "scenario", "method"))

#Estimate mean number of secondary cases divided by pi
out <- out %>% 
  mutate(R_sim_sup=round(R_sim/inf, 0),
         R_sim_inf=round(R_sim/sup, 0))

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
out$scen_sim <- case_when(out$scen_sim == "A1" ~ "CNTrW",
                          out$scen_sim == "B1" ~ "CTrW",
                          out$scen_sim == "S4" ~ "CTrW_H")

#Is the real contribution in the 95%HPD?
out <- out %>%
  mutate(correct=case_when(R<R_sim_inf ~ "Overestimation",
                           R>R_sim_sup ~ "Underestimation",
                           TRUE ~ "Correct"))

#------------------------------------
#Color code
strans <- c("Underestimation" = "steelblue",
            "Correct" = "chartreuse2",
            "Overestimation" = "indianred3")

#------------------------------------
#Select transmission scenario CTrW only
fig_1 <- out %>% filter(scen_sim=="CTrW") %>% 
  #filter(sp_infector=="cattle") #%>% 
  #filter(scenario!="SB" & scenario!="T+SB") 
  filter(scenario=="Reference")

fig_1$sp_infector <- factor(fig_1$sp_infector, levels=c("cattle", "badger", "boar"))

(A <- ggplot(data=fig_1)+
    geom_errorbar(aes(y=reorder_within(sim, R, within=sp_infector), xmin=R_sim_sup, xmax=R_sim_inf, color=correct))+
    geom_point(aes(y=reorder_within(sim, R, within=sp_infector), x= R, color = correct))+
    facet_grid(sp_infector~method, scales="free_y")+
    scale_x_log10()+
    scale_color_manual(name="Estimation", values=strans)+
    theme_bw(base_size=11)+
    labs(x="Number of transmission events", y=""))

A + theme(axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          legend.position="bottom")

ggsave("All_sp_ref_within.tiff",device="tiff", dpi="print", width=17, height=15, units="cm",  compression="lzw", scale=1)


#------------------------------------
#Select high mutation transmission scenario 
fig_2 <- out %>% filter(scen_sim=="CTrW_H") %>% filter(!is.na(sp_infector))

fig_2$sp_infector <- factor(fig_2$sp_infector, levels=c("cattle", "badger", "boar"))

(A <- ggplot(data=fig_2)+
    geom_errorbar(aes(y=reorder_within(sim, R, within=sp_infector), xmin=R_sim_sup, xmax=R_sim_inf, color=correct))+
    geom_point(aes(y=reorder_within(sim, R, within=sp_infector), x= R, color = correct))+
    facet_grid(sp_infector~method, scales="free_y")+
    scale_x_continuous()+
    scale_color_manual(name="Estimation", values=strans)+
    theme_bw(base_size=11)+
    labs(x="Number of transmission events", y=""))

A + theme(axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          legend.position="bottom")

ggsave("Host_contrib_high_mutation.tiff",device="tiff", dpi="print", width=17, height=15, units="cm",  compression="lzw", scale=1)

#------------------------------------
#Select dead-end transmission scenario 
fig_3 <- out %>% filter(scen_sim=="CNTrW") %>% filter(scenario=="Reference") %>%
  filter(!is.na(sp_infector))

fig_3$sp_infector <- factor(fig_3$sp_infector, levels=c("cattle", "badger", "boar"))

(A <- ggplot(data=fig_3)+
    geom_errorbar(aes(y=reorder_within(sim, R, within=sp_infector), xmin=R_sim_sup, xmax=R_sim_inf, color=correct))+
    geom_point(aes(y=reorder_within(sim, R, within=sp_infector), x= R, color = correct))+
    facet_grid(sp_infector~method, scales="free_y")+
    scale_x_continuous()+
    scale_color_manual(name="Estimation", values=strans)+
    theme_bw(base_size=11)+
    labs(x="Number of transmission events", y=""))

A + theme(axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          legend.position="bottom")

ggsave("Host_contrib_dead_end.tiff",device="tiff", dpi="print", width=17, height=15, units="cm",  compression="lzw", scale=1)
