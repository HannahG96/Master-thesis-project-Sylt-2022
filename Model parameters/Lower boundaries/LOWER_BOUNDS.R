############################### LOWER BOUNDARIES ########################################

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/Lower boundaries")

#load packages
library(writexl)
library(data.table)

seasons<-c("2009 spring","2009 summer","2009 autumn","2010 winter",
           "2010 spring","2010 summer","2010 autumn","2011 winter")

######################################################
### MINIMUM RESPIRATION OF CILIATES & MESOZOOPLANKTON
############################ Vézina & Platt, 1988  
minRESP<-as.data.frame(matrix(NA,nrow=2,ncol=(length(seasons)+1),dimnames=list(NULL,c("Comp",seasons))))
#Import seasonal SST estimates of SRB:
SST<-fread(file = "SST.csv", na.strings = "", dec = "," , data.table = FALSE)
minRESP$Comp<-c("Cil","Meso")
minRESP_Cil<-0.01 * 7.2 * exp(0.0693*SST[2:nrow(SST),"seasonal_AVG"])
minRESP_Meso<-0.01 * 2.3 * exp(0.0693*SST[2:nrow(SST),"seasonal_AVG"])
minRESP[1,2:ncol(minRESP)]<-minRESP_Cil
minRESP[2,2:ncol(minRESP)]<-minRESP_Meso

#Export minimum respiration rates as excel file:
#write_xlsx(minRESP, 
 #path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/minRESP.xlsx")


