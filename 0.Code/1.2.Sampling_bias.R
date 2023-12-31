rm(list=ls())
#------------------------------------
#package
library(tidyverse) #version 1.3.0

#------------------------------------
##Description:
##Create the set of sequences for each sampling scheme

#------------------------------------
#------------------------------------
#Transmission scenario considered
samp <- "B1" #B1 (reference scenario) or A1 (dead-end), B2 (badger index), S1 (single-host), S4 (high mutation rate)

for (j in 1:30){ #all trees in transmission scenario

  #Import sequence simulated from the reference tree
  seq_1 <- read_csv(paste0("./seq_",substr(samp, 1, 1),"/seq_",samp,"_",j,"_1.csv")) #no bias
  
  
  #------------------SAMPLING SCHEMES------------------------------------
  #seq_2: no wildlife before 2012 (T)
  seq_2 <- seq_1 %>% filter(sp == "boar" & tremoved >= 60 |
                        sp == "badger" & tremoved >= 60 |
                        sp == "cattle") #60 is the number of months from Jan 2007 until Jan 2012

  #seq_3: no boars (SW)
  seq_3 <- seq_1 %>% filter(sp == "badger" | sp == "cattle")
  

  #seq_4: no badgers (SB)
  seq_4 <- seq_1 %>% filter(sp == "boar" | sp == "cattle")
  

  #seq_5: no boars + no wildlife before 2012 (T+SW)
  seq_5 <- seq_2 %>% filter(sp == "badger" | sp == "cattle")

  #seq_6: no badgers + no wildlife before 2012 (T+SB)
  seq_6 <- seq_2 %>% filter(sp == "boar" | sp == "cattle")
  
  #output
  seq <- list(seq_1, seq_2, seq_3, seq_4, seq_5, seq_6)
  for (k in 2:6){
    write_csv(seq[[k]], paste0("seq_",samp,"_",j,"_",k,".csv"))
  }

  
}

