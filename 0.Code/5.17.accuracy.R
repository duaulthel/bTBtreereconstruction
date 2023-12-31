rm(list=ls())
#------------------------------------
#package
library(tidyverse) #version 1.3.0
library(lme4) #version 1.1-29
library(ggplot2) #version 3.3.5

#------------------------------------
##Description:
##Binomial GLMM on the accuracy indicator
##Code for Fig 2

#------------------------------------
#------------------------------------
#import accuracy results from every method
samp <- "B1" #B1 (reference scenario)

out_s1 <- read_csv(paste0("./Accuracy/seqTrack_",samp,"_acc.csv"))
out_s1$method <- "seqTrack"

out_o1 <- read_csv(paste0("./Accuracy/outbreaker_",samp,"_acc.csv"))
out_o1$method <- "outbreaker2"

out_t2 <- read_csv(paste0("./Accuracy/transphylo_WIW_",samp,"_acc.csv"))
out_t2$method <- "TransPhylo"

out <- rbind(out_s1, out_o1, out_t2)


#keep only trees that converged in TransPhylo
class(out$scenario) <- "character"
indice <- out %>% filter(method=="TransPhylo")
out <- out %>% filter(sim %in% indice$sim)

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
out$scen_sim <- "CTrW" #meaning cattle index and wild boars that transmit

#Random effect variable
out$sim <- as.factor(out$sim)

#------------------------------------
#Figure with results for Reference sampling scheme 
fig <- out %>% filter(scenario=="Reference") %>%
  mutate(prop=nb_acc/(nb_acc+nb_non_acc)*100)

(A <- ggplot(data=fig,aes(y=prop, x = method))+
    geom_boxplot()+
    theme_bw(base_size=11)+
    labs(x="Outbreak reconstruction method", y="Proportion of correctly identified
    transmission events (%)"))

ggsave("Accuracy_CW.tiff",device="tiff", dpi="print", width=17, height=10, units="cm",  compression="lzw", scale=1) #Fig 2

#------------------------------------
#Table with results from all methods (S2 Table)
tab <- out %>% 
  mutate(prop=nb_acc/(nb_acc+nb_non_acc)*100)

tab <- tab %>% group_by(scen_sim, scenario, method) %>%
  summarise(median=round(median(prop),1),
            min=round(min(prop),1),
            max=round(max(prop),1))

tab <- tab %>% mutate(res=paste0(median, " (", min, "-", max, ")")) %>%
  select(scen_sim, scenario, method, res)

tab <- tab %>% pivot_wider(names_from=scen_sim, values_from=res)

#------------------------------------
#Implement the GLMM on the accuracy results
model <- glmer(cbind(nb_acc, nb_non_acc) 
             ~ method + method:scenario + (1|sim),
             control = glmerControl(optimizer="bobyqa"),
             data=out,
             family=binomial)

summary(model) #and to get p-value

#Estimate 95% confidence interval
cc <- confint(model,parm="beta_",method="Wald")
ctab <- cbind(est=fixef(model),cc)

#Exponentiate results to get Odds Ratio
rtab <- exp(ctab)

#Create data.frame with results for Table 2
tab <- as.data.frame(round(rtab, 2))

