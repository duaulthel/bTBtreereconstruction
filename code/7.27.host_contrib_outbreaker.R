rm(list=ls())
#------------------------------------
#packages
library(tidyverse)

#------------------------------------
##Description:
##Estimate mean number of secondary cases 
##Comparison with the induced subtree

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

#If there is a sampling scheme different from A, we need the tree from induced_subtree
sup <- ifelse(substr(samp,1,1)=="S", "det", "det_sup")

#--------------------------------
for (j in 1:30){#j is the tree considered
  
  #initialize output
  nb <- NULL
  
  for (i in 1:nb_scheme){#sampling schemes
    
    sim <- paste0(samp,"_",j) 
    
    #reference tree
    Ttree <- read_csv(paste0("./Ttrees_",substr(ref_tree,1,1),"/Ttree_",sup,"_",ref_tree,"_",j,".csv")) 
    
    #if sampled scheme, take into account the induced subtree
    if(i>1){
      Ttree <- Ttree[(Ttree[[paste0("presence_",i)]])==1,]
    }
    
    #estimate number of hosts present in the induced subtree
    nb_hosts <- Ttree %>% group_by(sp_infected) %>%
      summarise(nb=n_distinct(infected)) %>% 
      rename(sp_infector=sp_infected)
    
    #remove index case (infector == NA) in order to only consider transmission events 
    Ttree <- Ttree %>% filter(!is.na(infector))
    
    #Estimate number of secondary cases per infector in the outbreak
    ref_nb <- Ttree %>% group_by(sp_infector) %>% 
      dplyr::count(infector) 
    
    #For each host-species sum the number of secondary cases
    ref_nb <- ref_nb %>% group_by(sp_infector) %>%
      summarise(R=sum(n))
    
    #Create new data.frame with number of secondary cases and number of hosts
    ref_nb <- full_join(ref_nb, nb_hosts, by=c("sp_infector"))
    
    ref_nb[is.na(ref_nb)] <- 0
    
    #--------------------------------
    #sampled sequences with host-species information
    seq <- read_csv(paste0("./seq_",substr(samp, 1, 1),"/seq_",sim,"_",i,".csv")) 
    seq$id <- 1:nrow(seq)#add id used in outbreaker2
    
    #results from outbreaker2
    tree <- read_csv(paste0("./outbreaker_",samp,"/res_tree",sim,"_",i,".csv")) 
    
    #merge tree with sequences with host-species information
    data <- merge(seq[,-4], tree, by.x="id", by.y="to")
    
    #find the name of infectors and their host-species
    data <- data %>% mutate(infector=infected[from],
                            sp_infector=sp[from])
    
    #Estimate number of sampled hosts
    nb_hosts_sim <- data %>% group_by(sp) %>%
      summarise(nb_sim=n_distinct(infected)) %>% 
      rename(sp_infector=sp)
    
    #remove index case (generations == NA) in order to consider direct transmission events
    data_sic <- data %>%
      filter(!is.na(generations)) 
    
    #Estimate number of secondary cases per infector in the reconstructed outbreak
    data_nb <- data_sic %>% group_by(sp_infector) %>% 
      dplyr::count(infector)
    
    #For each host-species sum the number of secondary cases
    data_nb <- data_nb %>% group_by(sp_infector) %>%
      summarise(R_sim=sum(n)) %>% filter(!is.na(sp_infector))
    
    #Create new data.frame with number of secondary cases and number of hosts
    data_nb <- full_join(data_nb, nb_hosts_sim, by=c("sp_infector"))
    
    data_nb[is.na(data_nb)] <- 0
    
    #--------------------------------
    #Create new data.frame with induced subtree and reconstructed tree results
    nb_trans <- full_join(ref_nb, data_nb, by=c("sp_infector"))
    
    nb_trans[is.na(nb_trans)] <- 0
    
    #Divide the sum of secondary cases by the number of hosts present
    nb_trans <- nb_trans%>% 
      mutate(R=R/nb, #for the induced subtree
             R_sim=R_sim/nb_sim, #for the reconstructed tree
             tree=sim, 
             scenario=i)
    
    nb <- rbind(nb, nb_trans)
  }
  
  #--------------------------------
  #Write output for transmission scenario
  if(j > 1){ #if this is not the first tree, call old file
    prev <- read_csv(paste0("outbreaker_",samp,"_reff.csv"))
    nb <- rbind(prev, nb)
    write_csv(nb, paste0("outbreaker_",samp,"_reff.csv"))
  }else{ #if this is the first tree, write new file
    write_csv(nb, paste0("outbreaker_",samp,"_reff.csv"))
  }
  
}
