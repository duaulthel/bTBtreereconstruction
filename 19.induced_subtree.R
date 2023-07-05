rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(igraph)

#------------------------------------
##Description:
##Estimate subtree induced from sampling scheme.

#------------------------------------
#------------------------------------
setwd("C:/Thèse_ANSES/codes_R/Ttrees_ref")

#Transmission scenario considered
samp <- "B1"


#------------------------------------

for (j in 1:30){#all trees in transmission scenario
  
  #reference tree
  Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_det_",ref_tree,"_",j,".csv")) 
  
  #Convert data.frame into igraph object
  tree_graph <- graph_from_data_frame(Ttree %>%
                                        filter(!is.na(infector)) %>%
                                        select(infector, infected), directed=TRUE)
  
  desc <- map(V(tree_graph), ~ names(subcomponent(tree_graph, .x, mode="in"))) %>%
    #subcomponent(tree_graph, v): breadth-first search is conducted starting from vertex v
    # Convert the list output to a data frame
    map_df(~data.frame(parents=.x), .id="child") %>% 
    filter(child != parents)
  
  #Add host-species and time of removal to the parents in the data.frame
  attr_desc <- left_join(desc, Ttree %>% 
                           select(infected, sp_infected, tremoved), 
                         by=c("parents"="infected")) %>% 
    rename(sp_infector=sp_infected, tremoved_inf=tremoved)
  attr_desc <- left_join(attr_desc, Ttree %>% 
                           select(infected, sp_infected, tremoved), 
                         by=c("child"="infected")) 
  
  #Find all infectors in the tree
  infectors <- unique(attr_desc$parents)
  
  #------------------------------------
  #Subtree induced by the temporal bias 
  subtree_T <- attr_desc %>%
    mutate(sampled=ifelse(sp_infected == "boar" & tremoved >= 60 |
                            sp_infected == "badger" & tremoved >= 60 |
                            sp_infected == "cattle", 1,0),
           sampled_inf=ifelse(sp_infector == "boar" & tremoved_inf >= 60 |
                                sp_infector == "badger" & tremoved_inf >= 60 |
                                sp_infector == "cattle", 1,0)) 
  
  #Remove unsampled infected hosts that did not transmit
  presence <- subtree_T %>% group_by(child) %>%
    summarise(presence_T=ifelse(sampled==0 & !(child %in% infectors), 0, 1))
  
  #Remove unsampled hosts that did not transmit to sampled hosts  
  presence_inf <- distinct(subtree_T %>% group_by(parents) %>%
                             summarise(presence_T_inf=ifelse(sampled_inf==0 & max(sampled)==0, 0, 1)) %>%
                             rename(child=parents))
  
  presence <- full_join(presence, presence_inf, by="child")
  presence[is.na(presence)] <- 1
  presence <- distinct(presence %>% 
                         mutate(presence_2=ifelse(presence_T_inf==0, 0, presence_T)))
  
  Ttree <- left_join(Ttree, presence %>% select(child, presence_2), by=c("infected"="child"))
  
  #------------------------------------
  #Subtree induced by the wild boar bias 
  subtree_W <- attr_desc %>% 
    mutate(sampled=ifelse(sp_infected == "badger" |
                            sp_infected == "cattle", 1, 0),
           sampled_inf=ifelse(sp_infector == "badger" |
                                sp_infector == "cattle", 1, 0))
  
  #Remove unsampled infected hosts that did not transmit
  presence <- subtree_W %>% group_by(child) %>%
    summarise(presence_W=ifelse(sampled==0 & !(child %in% infectors), 0, 1))
  
  #Remove unsampled hosts that did not transmit to sampled hosts  
  presence_inf <- distinct(subtree_W %>% group_by(parents) %>%
                             summarise(presence_W_inf=ifelse(sampled_inf==0 & max(sampled)==0, 0, 1)) %>%
                             rename(child=parents))
  
  presence <- full_join(presence, presence_inf, by="child")
  presence[is.na(presence)] <- 1
  presence <- distinct(presence %>% 
                         mutate(presence_3=ifelse(presence_W_inf==0, 0, presence_W)))
  
  Ttree <- left_join(Ttree, presence %>% select(child, presence_3), by=c("infected"="child"))
  
  #------------------------------------
  #Subtree induced by the badger bias 
  subtree_B <- attr_desc %>% 
    mutate(sampled=ifelse(sp_infected == "boar" |
                            sp_infected == "cattle", 1, 0),
           sampled_inf=ifelse(sp_infector == "boar" |
                                sp_infector == "cattle", 1, 0))
  #Remove unsampled infected hosts that did not transmit
  presence <- subtree_B %>% group_by(child) %>%
    summarise(presence_B=ifelse(sampled==0 & !(child %in% infectors), 0, 1))
  
  #Remove unsampled hosts that did not transmit to sampled hosts  
  presence_inf <- distinct(subtree_B %>% group_by(parents) %>%
                             summarise(presence_B_inf=ifelse(sampled_inf==0 & max(sampled)==0, 0, 1)) %>%
                             rename(child=parents))
  
  presence <- full_join(presence, presence_inf, by="child")
  presence[is.na(presence)] <- 1
  presence <- distinct(presence %>% 
                         mutate(presence_4=ifelse(presence_B_inf==0, 0, presence_B)))
  
  Ttree <- left_join(Ttree, presence %>% select(child, presence_4), by=c("infected"="child"))
  
  #------------------------------------
  #Subtree induced by the temporal + Wild boar bias 
  subtree_TW <- attr_desc %>%
    mutate(sampled=ifelse(sp_infected == "badger" & tremoved >= 60 |
                            sp_infected == "cattle", 1,0),
           sampled_inf=ifelse(sp_infector == "badger" & tremoved_inf >= 60 |
                                sp_infector == "cattle", 1,0)) 
  
  #Remove unsampled infected hosts that did not transmit
  presence <- subtree_TW %>% group_by(child) %>%
    summarise(presence_TW=ifelse(sampled==0 & !(child %in% infectors), 0, 1))
  
  #Remove unsampled hosts that did not transmit to sampled hosts  
  presence_inf <- distinct(subtree_TW %>% group_by(parents) %>%
                             summarise(presence_TW_inf=ifelse(sampled_inf==0 & max(sampled)==0, 0, 1)) %>%
                             rename(child=parents))
  
  presence <- full_join(presence, presence_inf, by="child")
  presence[is.na(presence)] <- 1
  presence <- distinct(presence %>% 
                         mutate(presence_5=ifelse(presence_TW_inf==0, 0, presence_TW)))
  
  Ttree <- left_join(Ttree, presence %>% select(child, presence_5), by=c("infected"="child"))
  
  #------------------------------------
  #Subtree induced by the temporal + Badger bias 
  subtree_TB <- attr_desc %>%
    mutate(sampled=ifelse(sp_infected == "boar" & tremoved >= 60 |
                            sp_infected == "cattle", 1,0),
           sampled_inf=ifelse(sp_infector == "boar" & tremoved_inf >= 60 |
                                sp_infector == "cattle", 1,0)) 
  
  #Remove unsampled infected hosts that did not transmit
  presence <- subtree_TB %>% group_by(child) %>%
    summarise(presence_TB=ifelse(sampled==0 & !(child %in% infectors), 0, 1))
  
  #Remove unsampled hosts that did not transmit to sampled hosts  
  presence_inf <- distinct(subtree_TB %>% group_by(parents) %>%
                             summarise(presence_TB_inf=ifelse(sampled_inf==0 & max(sampled)==0, 0, 1)) %>%
                             rename(child=parents))
  
  presence <- full_join(presence, presence_inf, by="child")
  presence[is.na(presence)] <- 1
  presence <- distinct(presence %>% 
                         mutate(presence_6=ifelse(presence_TB_inf==0, 0, presence_TB)))
  
  Ttree <- left_join(Ttree, presence %>% select(child, presence_6), by=c("infected"="child"))
  
  #------------------------------------
  #Write output for detailed tree with presence (1/0) in the induced subtree
  write_csv(Ttree, paste0("Ttree_det_sup_",samp,"_",j,".csv"))
}  

