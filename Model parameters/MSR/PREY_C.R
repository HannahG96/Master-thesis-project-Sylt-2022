#PREY CONCENTRATIONS TO CALCULATE JELLY FISH INGESTION RATES
#--> IR=CR*C_PREY
#-->use of median Noctiluca & diatom biomass to limit prey concentration and get realistic ingestion estimates

#set working directory:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/MSR")

#load packages:
library(writexl)
library(data.table)
#Define prey species of jelly fish:
prey.Hydro<-c("Cil","Tun","Clado","Cop","Biv","Gastr")
prey.Ppil<-c("Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro")
prey.Mlei<-c("Dia","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro")
prey.Bcu<-c("Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Ppil","Mlei")

#Import seasonal carbon biomass of prey species:
PREY_C<-fread(file = "CBIOM.csv", na.strings = "", dec = "," , data.table = FALSE)
C_Nsci<-fread(file = "CBIOM_Nsci.csv", na.strings = "", dec = "," , data.table = FALSE)#seasonal median biomass of Noctiluca
C_Nsci<-C_Nsci[-which(C_Nsci$season=="2009 winter"),]
C_Dia<-fread(file = "CBIOM_Dia.csv", na.strings = "", dec = "," , data.table = FALSE)#seasonal median biomass of diatoms
PREY_C<-merge(PREY_C,C_Nsci[,c("season","seasonal_MED")],by="season")
colnames(PREY_C)[c(4,19)]<-c("Nsci.avg","Nsci")
PREY_C<-merge(PREY_C,C_Dia[,c("season","seasonal_MED")],by="season")
colnames(PREY_C)[c(3,20)]<-c("Dia.avg","Dia")

#CHECK:
#calculate prey concentration using Nsci & Dia median biomass:
prey_C.Hydro<-mean(apply(PREY_C[,prey.Hydro],1,sum,na.rm=TRUE))#12.80317
prey_C.Mlei<-mean(apply(PREY_C[,prey.Mlei],1,sum,na.rm=TRUE))#96.53422
prey_C.Ppil<-mean(apply(PREY_C[,prey.Ppil],1,sum,na.rm=TRUE))#41.42698
prey_C.Bcu<-mean(apply(PREY_C[,prey.Bcu],1,sum,na.rm=TRUE))#41.36948

#Store prey concentrations of jelly fish in data frame
seasons<-c(paste(2009,c("spring","summer","autumn"),sep=" "),
           paste(2010,c("winter","spring","summer","autumn"),sep=" "),"2011 winter")
preys<-list(prey.Hydro,prey.Ppil,prey.Mlei,prey.Bcu)
PREY_C.jellies<-as.data.frame(matrix(NA,nrow=4,ncol=(length(seasons)+1),
                                     dimnames=list(NULL,c("Compartment",seasons))))
PREY_C.jellies$Compartment<-c("Hydro","Ppil","Mlei","Bcu")
for(i in 1:length(preys)){
  for(u in 2:ncol(PREY_C.jellies)){
  row<-which(PREY_C$season==seasons[u-1])
  PREY_C.jellies[i,u]<-apply(PREY_C[row,preys[[i]]],1,sum)}}
PREY_C.jellies$MEAN<-c(prey_C.Hydro,prey_C.Ppil,prey_C.Mlei,prey_C.Bcu)

#Export jellyfish prey biomass concentrations as Excel file:
#write_xlsx(PREY_C.jellies, 
          #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/PREY_C.xlsx")
           


