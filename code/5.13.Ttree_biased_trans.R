rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(igraph)

#------------------------------------
##Description:
##Add columns (one per bias scenario) with transmission events that 
##take into account missing hosts in the reference trees 

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1"


for (j in 1:30){#all trees in the transmission scenario
  
  #reference tree
  Ttree <- read_csv(paste0("./Ttrees_",substr(samp,1,1),"/Ttree_",samp,"_",j,".csv")) 
  
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
  attr_desc <- attr_desc %>% mutate(trans = paste(parents, child, sep="->"))
  
  #------------------------------------
  #Sampling scheme A
  trans_1 <- attr_desc %>% group_by(child) %>% 
    slice_head(n=1) %>% #select closest parent
    select(child, trans) %>% rename(trans_1=trans)
  
  Ttree <- left_join(Ttree, trans_1, by=c("infected"="child"))
  
  #------------------------------------
  #Sampling scheme T
  trans_2 <- attr_desc %>% group_by(child) %>%
    filter(sp_infected == "boar" & tremoved >= 60 |
             sp_infected == "badger" & tremoved >= 60 |
             sp_infected == "cattle") %>%
    filter(sp_infector == "boar" & tremoved_inf >= 60 |
             sp_infector == "badger" & tremoved_inf >= 60 |
             sp_infector == "cattle") %>% slice_head(n=1) %>%
    select(child, trans)%>% rename(trans_2=trans)
  
  Ttree <- left_join(Ttree, trans_2, by=c("infected"="child"))
  
  #------------------------------------
  #Sampling scheme SW
  trans_3 <- attr_desc %>% group_by(child) %>%
    filter(sp_infected == "badger" |
             sp_infected == "cattle") %>%
    filter(sp_infector == "badger" |
             sp_infector == "cattle")  %>% slice_head(n=1) %>%
    select(child, trans)%>% rename(trans_3=trans)
  
  Ttree <- left_join(Ttree, trans_3, by=c("infected"="child"))
  
  #------------------------------------
  #Sampling scheme SB
  trans_4 <- attr_desc %>% group_by(child) %>%
    filter(sp_infected == "boar" |
             sp_infected == "cattle") %>%
    filter(sp_infector == "boar" |
             sp_infector == "cattle")  %>% slice_head(n=1) %>%
    select(child, trans)%>% rename(trans_4=trans)
  
  Ttree <- left_join(Ttree, trans_4, by=c("infected"="child"))
  
  #------------------------------------
  #Sampling scheme T+SW
  trans_5 <- attr_desc %>% group_by(child) %>%
    filter(sp_infected == "badger" & tremoved >= 60 |
             sp_infected == "cattle") %>%
    filter(sp_infector == "badger" & tremoved_inf >= 60 |
             sp_infector == "cattle") %>% slice_head(n=1) %>%
    select(child, trans)%>% rename(trans_5=trans)
  
  Ttree <- left_join(Ttree, trans_5, by=c("infected"="child"))
  
  #------------------------------------
  #Sampling scheme T+SB
  trans_6 <- attr_desc %>% group_by(child) %>%
    filter(sp_infected == "boar" & tremoved >= 60 |
             sp_infected == "cattle") %>%
    filter(sp_infector == "boar" & tremoved_inf >= 60 |
             sp_infector == "cattle") %>% slice_head(n=1) %>%
    select(child, trans)%>% rename(trans_6=trans)
  
  Ttree <- left_join(Ttree, trans_6, by=c("infected"="child"))
  
  #------------------------------------
  #Write output
  write_csv(Ttree, paste0("Ttree_det_",samp,"_",j,".csv"))
}  
