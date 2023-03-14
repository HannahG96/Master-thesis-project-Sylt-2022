################################# SPECIFIC INGESTION ##############################

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/")

#load packages
library(varhandle)
library(writexl)

seasons<-c("2009 spring","2009 summer","2009 autumn","2010 winter",
           "2010 spring","2010 summer","2010 autumn","2011 winter")

#########################################################
### HERRING DIET {fraction of total ingestion}
######################################## Kellnreitner et al. (2013)
#-->Use dietary range for all seasons (???)
DIET_Her<-data.frame(season=c("MIN summer", "MAX summer", "autumn"),Dia=NA,Tun=NA,Clado=NA,
                     Cop=NA,Biv=NA,Gastr=NA,Poly=NA,Hydro=NA,Other=NA)
DIET_Her$Dia<-c(0,0,0.017)
DIET_Her$Tun<-c(0,0.084,0)
DIET_Her$Clado<-c(0,0.004,0)
DIET_Her$Cop<-c((0.256+0.122),(0.666+0.183),(0.518+0.072))
DIET_Her$Biv<-c(0.001,0.013,0.001)
DIET_Her$Gastr<-c(0.006,0.081,0)
DIET_Her$Poly<-c(0.008,0.053,0)
DIET_Her$Hydro<-c(0,0.02,0)
0.064+0.306+0.007+0.092 #June
0.009+0.03+0.031+0.012+0.018 #July
0.058+0.056+0.045+0.064+0.013 #Aug
0.081+0.109+0.1+0.099 #Sept
DIET_Her$Other<-c(0.1,0.469,0.389)
#Export Herring diet as excel file:
#write_xlsx(DIET_Her, 
    # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/DIET_Her.xlsx")

#########################################################
### M.LEIDYI DIET {fraction of total ingestion}
######################################## Kellnreitner et al. (2013)
#-->Use dietary range for all seasons (???)
DIET_Mlei<-data.frame(season=c("MIN summer", "MAX summer", "autumn"),Dia=NA,Tun=NA,Nsci=NA,
                      Clado=NA,Cop=NA,Biv=NA,Gastr=NA,Poly=NA,Hydro=NA,Other=NA)
DIET_Mlei$Dia<-c(0,0,0.066)
DIET_Mlei$Nsci<-c(0,0.35,0.012)
DIET_Mlei$Tun<-c(0,0.042,0.017)
DIET_Mlei$Clado<-c(0,0.068,0)
DIET_Mlei$Cop<-c((0.111+0),(0.497+0.125),(0.353+0.072))
DIET_Mlei$Biv<-c(0,0.472,0.155)
DIET_Mlei$Gastr<-c(0,0.029,0.017)
DIET_Mlei$Poly<-c(0,0.048,0)
DIET_Mlei$Hydro<-c(0,0.097,0)
0.0038+0.213+0.003 #June
0.036+0.089+0.004+0.041+0.036+0.037+0.03 #July
0.03+0.107+0.034+0.049+0.031 #Aug
0.047+0.025+0.164+0.064#Sept
DIET_Mlei$Other<-c(0.2198,0.273,0.3)

#Export M.leidyi diet as excel file:
#write_xlsx(DIET_Mlei, 
 #path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/DIET_Mlei.xlsx")

##############################################################
### SPECIFIC CLEARANCE OF DIATOMS/MICROZOOPLANKTON BY COPEPODS {cm3/idv/hr}
######################################## Gasparini (2000)
#1.Define specific Weight of copepod species: A. clausii, C. hamatus, T.longicornis
W_A.clausii<-0.001507118 #{mg C/idv}
W_C.hamatus<-0.004359243 #{mg C/idv}
W_T.longicornis<-0.008740312 #{mg C/idv}

#2.Get microplankton and copepod carbon contributions
source("Phyto-_Microplankton.spp_J.Rick/format_PhytoMicroplanktonData.R")
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/")
source("Zooplankton_P.Martens/format_ZooplanktonData.R")
rm.all.but(keep=c("seasons","TOTCARB.SRB","CARBcontr.ZOO","W_A.clausii","W_C.hamatus",
                  "W_T.longicornis"), 
           envir=.GlobalEnv, keep_functions=FALSE, gc_limit=100,regex="auto")

#Format microplankton data
#add season column:
TOTCARB.SRB$yearS<-TOTCARB.SRB$year
TOTCARB.SRB$season<-NA
TOTCARB.SRB[which(TOTCARB.SRB$month==12),"yearS"]<-TOTCARB.SRB[which(TOTCARB.SRB$month==12),"yearS"]+1
for(i in 1:nrow(TOTCARB.SRB)){
  if(length(which(c(12,1,2)==TOTCARB.SRB[i,"month"]))>0)TOTCARB.SRB[i,"season"]<-"winter"
  if(length(which(c(3,4,5)==TOTCARB.SRB[i,"month"]))>0)TOTCARB.SRB[i,"season"]<-"spring"
  if(length(which(c(6,7,8)==TOTCARB.SRB[i,"month"]))>0)TOTCARB.SRB[i,"season"]<-"summer"
  if(length(which(c(9,10,11)==TOTCARB.SRB[i,"month"]))>0)TOTCARB.SRB[i,"season"]<-"autumn"}
TOTCARB.SRB$season<-paste(TOTCARB.SRB$yearS,TOTCARB.SRB$season,sep=" ")
#Select time period spring 2009-winter 2011 and re-order rows chronologically
for(i in 1:length(seasons)){
  if(i==1){
    TOTCARB<-TOTCARB.SRB[which(TOTCARB.SRB$season==seasons[i]),]
  }else{
    TOTCARB<-rbind(TOTCARB,TOTCARB.SRB[which(TOTCARB.SRB$season==seasons[i]),])}}
#Calculate seasonal diatom (only those included in model) + microzooplankton 
#(total biomass, incl. ciliates, het. flagellates, dinoflagellates....) concentrations:
#in {mg C/m3}
PREY_C<-data.frame(season=seasons,Diatoms=NA,Microzooplankton=NA)
for(i in 1:length(seasons)){
  mydat<-TOTCARB[which(TOTCARB[,"season"]==seasons[i]),]
  PREY_C[i,"Diatoms"]<-mean(unname(tapply(mydat$Diatoms,INDEX=mydat$month,FUN=mean)))
  PREY_C[i,"Microzooplankton"]<-mean(unname(tapply(mydat$Protozooplankton.SRB,
                                                   INDEX=mydat$month,FUN=mean)))}
    
#3.Calculate seasonal weight-specific clearancec rate {m3/mg C/day} of phyto-/microzooplankton by copepods:
#CR from Gasparini (2000) in cm3/idv/hr
CR_Cop<-data.frame(Cop=c("Acartia_sp","C.hamatus","T.longicornis","AVG"), 
                   minCR_Dia=c(0.6,0.5,0.8,mean(0.6,0.5,0.8)),
                   maxCR_Dia=c(1.1,0.8,0.8,mean(1.1,0.8,0.8)),
                   minCR_Mic=c(0.5,1.2,0.3,mean(0.5,1.2,0.3)),
                   maxCR_Mic=c(1.4,1.8,1,mean(1.4,1.8,1)))
CR_Cop$minDia_m3mgCd<-NA
CR_Cop$maxDia_m3mgCd<-NA
CR_Cop$minMic_m3mgCd<-NA
CR_Cop$maxMic_m3mgCd<-NA
#CR from Gasparini (2000) in m3/mg C/day
CR_Cop[1,6:9]<-CR_Cop[1,2:5]/W_A.clausii*0.000001*12
CR_Cop[2,6:9]<-CR_Cop[2,2:5]/W_C.hamatus*0.000001*12
CR_Cop[3,6:9]<-CR_Cop[3,2:5]/W_T.longicornis*0.000001*12
CR_Cop[4,6:9]<-CR_Cop[4,2:5]/mean(W_A.clausii,W_C.hamatus,W_T.longicornis)*0.000001*12

CARBcontr_A.clausii<-CARBcontr.ZOO[which(CARBcontr.ZOO[,"Species"]=="Acartia_sp"),]
CARBcontr_C.hamatus<-CARBcontr.ZOO[which(CARBcontr.ZOO[,"Species"]=="C.hamatus"),]
CARBcontr_T.longicornis<-CARBcontr.ZOO[which(CARBcontr.ZOO[,"Species"]=="T.longicornis"),]

SpecING_Cop<-data.frame(season=seasons,PREY_C.Dia=PREY_C$Diatoms,minCR_Dia=NA,maxCR_Dia=NA,
                        PREY_C.Mic=PREY_C$Microzooplankton,minCR_Mic=NA,maxCR_Mic=NA)
for(i in 1:length(seasons)){
  CARBcontr<-rep(NA,4)
  CARBcontr[1]<-CARBcontr_A.clausii[which(CARBcontr_A.clausii[,"Season"]==seasons[i]),"pct_Carbon"]/100
  CARBcontr[2]<-CARBcontr_C.hamatus[which(CARBcontr_C.hamatus[,"Season"]==seasons[i]),"pct_Carbon"]/100
  CARBcontr[3]<-CARBcontr_T.longicornis[which(CARBcontr_T.longicornis[,"Season"]==seasons[i]),"pct_Carbon"]/100
  CARBcontr[4]<-1-CARBcontr[1]-CARBcontr[2]-CARBcontr[3]
  SpecING_Cop[i,"minCR_Dia"]<-sum(CR_Cop$minDia_m3mgCd*CARBcontr)
  SpecING_Cop[i,"maxCR_Dia"]<-sum(CR_Cop$maxDia_m3mgCd*CARBcontr)
  SpecING_Cop[i,"minCR_Mic"]<-sum(CR_Cop$minMic_m3mgCd*CARBcontr)
  SpecING_Cop[i,"maxCR_Mic"]<-sum(CR_Cop$maxMic_m3mgCd*CARBcontr)}

#Export specific ingestion rates of copepods as Excel file:
#write_xlsx(SpecING_Cop, 
          #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/SpecING_Cop.xlsx")

######################################################################
### SPECIFIC INGESTION BY BIVALVE LARVAE {fraction of total ingestion}
############################################## Lindeque et al. (2015)  
SpecING_Biv<-data.frame(prey=c("Dia+Phae+Imp","Phae","Cil","Bac"),
                        MIN=c(0.75,0,0.006,NA),MAX=c(NA,0.7,0.025,0.02))
#Export specific ingestion rates of bivalves as Excel file:
write_xlsx(SpecING_Biv, 
 path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/SpecING_Biv.xlsx")
