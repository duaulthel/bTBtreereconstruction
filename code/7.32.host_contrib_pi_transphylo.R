rm(list=ls())
#------------------------------------
#packages
library(tidyverse)
library(lubridate)

#------------------------------------
##Description:
##Estimate number of transmission events due to each host-species 
##Comparison with the reference tree

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1"

#------------------------------------
non_se <- c("A1", "B1", "B2", "S1") 
#tree_name of reference tree same as reconstructed tree name

#For the reconstructed trees with a higher mutation rate
#Find name of the reference tree
ref_tree <- case_when(samp %in% non_se ~ samp,
                      samp=="S4" ~ "B1")

nb_scheme <- ifelse(samp=="B1", 6, 1)

#--------------------------------
#Not all trees converged, list of those that did:
# samp=="A1" ~ c(1,6,7,8,10,13,16,17,20,21,22,24,25,26,27,28,29)
# samp=="B1" ~ c(1:4,8,10:14,17,18,20:23,26:30)
# samp=="B2" ~ c(1:6,9,15:18,20,22:26,29)
# samp=="S1" ~ c(1:2,4,6:15,17:26,28:30)
# samp=="S4" ~ c(1:3,6,8:9,11:15,17:18,20:22,24:27,29)

#----------------------------------------------------------------------------

for (j in c(1:2,4:6,8:9,11:16,18:19,22:29)){#j is the tree considered
  sim <- paste0(samp,"_",j) 
  
  Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_det_",ref_tree,"_",j,".csv")) #reference tree
  
  #estimate number of hosts present
  nb_hosts <- Ttree %>% group_by(sp_infected) %>%
    summarise(nb=n_distinct(infected)) %>% 
    rename(sp_infector=sp_infected)
  
  #remove index case (infector == NA) in order to only consider transmission events 
  Ttree <- Ttree %>% filter(!is.na(infector))
  
  #Estimate number of secondary cases per infector in the outbreak
  ref_nb <- Ttree %>% group_by(sp_infector)  %>% 
    dplyr::count(infector)
  
  #For each host-species sum the number of secondary cases
  ref_nb <- ref_nb %>% group_by(sp_infector) %>%
    summarise(R=sum(n))
  
  #Create new data.frame with number of secondary cases and number of hosts
  ref_nb <- full_join(ref_nb, nb_hosts, by=c("sp_infector"))
  
  ref_nb[is.na(ref_nb)] <- 0
  
  #--------------------------------
  #initialize output
  nb <- NULL
  
  for (i in 1:nb_scheme){#sampling schemes
    #sampled sequences with host-species information
    seq <- read_csv(paste0("./seq_",substr(samp, 1, 1),"/seq_",sim,"_",i,".csv"))
    tstart <- as.Date("01/01/07", "%d/%m/%y")
    seq$tremoved <- tstart %m+% months(seq$tremoved)
    seq <- seq %>% 
      mutate(infected_id = paste(infected,sp,tremoved, sep="_")) 
    #in order to have the same name as in beast
    
    #results from medoid TransPhylo
    tree <- read_csv(paste0("./TransPhylo_",samp,"/ttree_",sim,"_",i,".csv"))
    
    data <- full_join(seq[,c(1,3,5)], tree, by="infected_id")
    
    #find the name of infectors and their host-species
    data <- data %>% mutate(infector=ifelse(infector_id!=0,infected[infector_id], "NA"),
                            sp_infector=ifelse(infector_id!=0,sp[infector_id], "NA"))
    
    #Estimate number of sampled hosts
    nb_hosts_sim <- data %>% group_by(sp) %>%
      summarise(nb_sim=n_distinct(infected)) %>% 
      rename(sp_infector=sp) %>% filter(!is.na(sp_infector))
    
    #remove index case and unsampled cases in order to consider direct transmission events
    data_sic <- data %>%
      filter(infector_id!=0) %>%
      filter(!(is.na(sp)))
    
    #Estimate number of secondary cases per infector in the reconstructed outbreak
    data_nb <- data_sic %>% group_by(sp_infector) %>% 
      dplyr::count(infector)
    
    #For each host-species sum the number of secondary cases
    data_nb <- data_nb %>% group_by(sp_infector) %>%
      summarise(R_sim=sum(n)) %>% filter(!(is.na(sp_infector)))
    
    #Create new data.frame with number of secondary cases and number of hosts
    data_nb <- full_join(data_nb, nb_hosts_sim, by=c("sp_infector"))
    
    data_nb[is.na(data_nb)] <- 0
    
    #--------------------------------
    #Create new data.frame with induced subtree and reconstructed tree results
    nb_trans <- full_join(ref_nb, data_nb, by=c("sp_infector"))
    
    nb_trans[is.na(nb_trans)] <- 0
    
    #Divide the sum of secondary cases by the number of hosts present
    nb_trans <- nb_trans%>% 
      mutate(R=R/nb, #for the reference tree
             R_sim=ifelse(nb_sim!=0, R_sim/nb_sim, 0), #for the reconstructed tree
             tree=sim, 
             scenario=i)
    
    nb <- rbind(nb, nb_trans)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("transphylo_",samp,"_reff_pi.csv"))
    nb <- rbind(prev, nb)
    write_csv(nb, paste0("transphylo_",samp,"_reff_pi.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(nb, paste0("transphylo_",samp,"_reff_pi.csv"))
  }
  
}
