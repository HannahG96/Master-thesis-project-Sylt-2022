####################### CHECK ENERGY BUDGETS OF COMPARTMENTS #############################

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/Parameters/")

#load packages
library(varhandle)
library(writexl)
library(data.table)

seasons<-c("2009 spring","2009 summer","2009 autumn","2010 winter",
           "2010 spring","2010 summer","2010 autumn","2011 winter")

#Import parameter files:
parameters<-c("BacGrowthEff","CBIOM","DIET_Her","DIET_Mlei","DOCexs","DOCuptake","GPP_NPP",
             "GROWTH_Phae","MCR","MGR","MIR","MRR","PhyGRAZING","POC_EXPORTtoNorthSea",
             "PREY_C","PtoR","SpecING_Cop","SpecING_Biv","SST")
PARAMS<-list()
files<-paste(parameters,".csv",sep="")
for(i in 1:length(files)){
  PARAMS[[i]]<-fread(file = files[i], na.strings = "", dec = "," , data.table = FALSE)}
names(PARAMS)<-parameters

#####################################
####PHYTOPLANKTON
#####################################
#-->Nutrient uptake
#-->Production of POC and DOC
#-->Phaeocystis production (growth)
#-->DOC exsudation
#-->Respiration
#-->Grazing by Microzooplankton
#-->Mortality
paramsPHY<-c("BIOM_Dia","BIOM_Phae","NutrUp","PROD_PocDoc","PROD_Phae","DOCexs","RESP","MicGRAZ",
             "MORT")
MIN_Phy<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsPHY)+1),
                                       dimnames=list(NULL,c("season",paramsPHY))))
MIN_Phy$season<-seasons
MAX_Phy<-MIN_Phy
MAX_Phy[,c("BIOM_Dia","BIOM_Phae")]<-MIN_Phy[,c("BIOM_Dia","BIOM_Phae")]<-PARAMS[["CBIOM"]][,c("Dia","Phae")]
MAX_Phy[,c("NutrUp","PROD_PocDoc","RESP")]<-PARAMS[["GPP_NPP"]][,c("maxGPP","maxNPP","maxGPPminusNPP")]
MIN_Phy[,c("NutrUp","PROD_PocDoc","RESP")]<-PARAMS[["GPP_NPP"]][,c("minGPP","minNPP","minGPPminusNPP")]
MAX_Phy[,"PROD_Phae"]<-PARAMS[["GROWTH_Phae"]][,"MAX {per day}"]*MAX_Phy$BIOM_Phae
MIN_Phy[,"PROD_Phae"]<-PARAMS[["GROWTH_Phae"]][,"MIN {per day}"]*MIN_Phy$BIOM_Phae
MAX_Phy[,"DOCexs"]<-PARAMS[["DOCexs"]][,"maxDOCexs"]*MAX_Phy$PROD_PocDoc
MIN_Phy[,"DOCexs"]<-PARAMS[["DOCexs"]][,"minDOCexs"]*MIN_Phy$PROD_PocDoc
MAX_Phy[,"MicGRAZ"]<-PARAMS[["PhyGRAZING"]][,"maxGRAZING"]*(MAX_Phy$BIOM_Dia+MAX_Phy$BIOM_Phae)
MIN_Phy[,"MicGRAZ"]<-PARAMS[["PhyGRAZING"]][,"minGRAZING"]*(MIN_Phy$BIOM_Dia+MIN_Phy$BIOM_Phae)
#Mortality:
#Between 5 and 50% of NPP (Vézina & Platt, 1988) 
MAX_Phy$MORT<-0.5*MAX_Phy$PROD_PocDoc
MIN_Phy$MORT<-0.05*MAX_Phy$PROD_PocDoc
#What production is available for grazing (part+dissPROD-R-DOCexs-MORT)?
MAX_Phy$PROD<-MAX_Phy$NutrUp-MAX_Phy$MORT-MAX_Phy$RESP-MAX_Phy$DOCexs
MIN_Phy$PROD<-MIN_Phy$NutrUp-MIN_Phy$MORT-MIN_Phy$RESP-MIN_Phy$DOCexs
#What production is left for diatoms?
MAX_Phy$PROD-MAX_Phy$PROD_Phae
MIN_Phy$PROD-MIN_Phy$PROD_Phae
#NOTES:
#ACCEPTED
#although VERY HIGH DOC exsudation flow: maybe DOC exsudation should be defined as fraction of PARTICULATE production

#####################################
####BACTERIA
#####################################
#-->DOC uptake
#-->Production
#-->Respiration
#-->Mortality
paramsBac<-c("DocUp","PROD","RESP","MORT")
MIN_Bac<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsBac)+1),
                              dimnames=list(NULL,c("season",paramsBac))))
MIN_Bac$season<-seasons
MAX_Bac<-MIN_Bac
#DOC uptake at maximum DOC exsudation:
MAX_Bac[,"DocUp"]<-PARAMS[["DOCuptake"]][,"DOCuptake"]*MAX_Phy$DOCexs
#DOC uptake at minimum DOC exsudation:
MIN_Bac[,"DocUp"]<-PARAMS[["DOCuptake"]][,"DOCuptake"]*MIN_Phy$DOCexs
#Respiration according to P/R and average phytoplankton PROD:
MAX_Bac[,"RESP"]<-apply(data.frame(min=MIN_Phy$PROD,max=MAX_Phy$PROD),1,mean)/PARAMS[["PtoR"]][,"maxPtoR"]
MIN_Bac[,"RESP"]<-apply(data.frame(min=MIN_Phy$PROD,max=MAX_Phy$PROD),1,mean)/PARAMS[["PtoR"]][,"minPtoR"]
#Mortality according to Vézina & Savenkoff (1999) and average Doc Uptake+respiration:
MAX_Bac[,"MORT"]<-0.75*apply(data.frame(max=MAX_Bac$DocUp,min=MIN_Bac$DocUp),1,mean)-
                               apply(data.frame(max=MAX_Bac$RESP,min=MIN_Bac$RESP),1,mean)
MIN_Bac[,"MORT"]<-0.5*apply(data.frame(max=MAX_Bac$DocUp,min=MIN_Bac$DocUp),1,mean)-
                               apply(data.frame(max=MAX_Bac$RESP,min=MIN_Bac$RESP),1,mean)
#Production according to BGE and assumed mean DOC uptake 
#(use BGE=BP/(DocUptake-MORT) instead of BGE=BP/(BP+R)):
MAX_Bac[,"PROD"]<-PARAMS[["BacGrowthEff"]][,"maxBGE"]*
  (apply(data.frame(max=MAX_Bac$DocUp,min=MIN_Bac$DocUp),1,mean)-
          apply(data.frame(max=MAX_Bac$MORT,min=MIN_Bac$MORT),1,mean))
MIN_Bac[,"PROD"]<-PARAMS[["BacGrowthEff"]][,"minBGE"]*
  (apply(data.frame(max=MAX_Bac$DocUp,min=MIN_Bac$DocUp),1,mean)-
     apply(data.frame(max=MAX_Bac$MORT,min=MIN_Bac$MORT),1,mean))
#Mass balance?
MAX_Bac$DocUp-MAX_Bac$PROD-MAX_Bac$RESP-MAX_Bac$MORT
MIN_Bac$DocUp-MIN_Bac$PROD-MIN_Bac$RESP-MIN_Bac$MORT
#TO BE FIGURED OUT IN LIM!

#####################################
####CILIATES
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsCil<-c("BIOM_Cil","PROD","ING","RESP")
MIN_Cil<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsCil)+1),
                              dimnames=list(NULL,c("season",paramsCil))))
MIN_Cil$season<-seasons
MAX_Cil<-MIN_Cil
MAX_Cil[,"BIOM_Cil"]<-MIN_Cil[,"BIOM_Cil"]<-PARAMS[["CBIOM"]][,"Cil"]
MAX_Cil[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Cil"),
                                  2:ncol(PARAMS[["MGR"]])])))*MAX_Cil$BIOM_Cil
MAX_Cil[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Cil"),
                                  2:ncol(PARAMS[["MIR"]])])))*MAX_Cil$BIOM_Cil
MAX_Cil[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Cil"),
                                                     2:ncol(PARAMS[["MRR"]])])))*MAX_Cil$BIOM_Cil
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Cil$RESP<-0.01*7.2 * exp(0.0693*PARAMS[["SST"]][2:nrow(PARAMS[["SST"]]),"seasonal_AVG"]) * MIN_Cil$BIOM_Cil
MAX_Cil$RESP-MIN_Cil$RESP 
#Limits of egestion:
MIN_Cil$EGEST<-0.1*MAX_Cil$ING #checked for highest possible ingestion
MAX_Cil$EGEST<-MAX_Cil$RESP
MAX_Cil$EGEST-MIN_Cil$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Cil$ING-MIN_Cil$RESP-MIN_Cil$EGEST
#NOTES:
#ACCEPTED

#####################################
####NOCTILUCA
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsNsci<-c("BIOM_Nsci","PROD","ING","RESP")
MIN_Nsci<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsNsci)+1),
                              dimnames=list(NULL,c("season",paramsNsci))))
MIN_Nsci$season<-seasons
MAX_Nsci<-MIN_Nsci
MAX_Nsci[,"BIOM_Nsci"]<-MIN_Nsci[,"BIOM_Nsci"]<-PARAMS[["CBIOM"]][,"Nsci"]
MAX_Nsci[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Nsci"),
                                                     2:ncol(PARAMS[["MGR"]])])))*MAX_Nsci$BIOM_Nsci
MAX_Nsci[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Nsci"),
                                                    2:ncol(PARAMS[["MIR"]])])))*MAX_Nsci$BIOM_Nsci
MAX_Nsci[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Nsci"),
                                                     2:ncol(PARAMS[["MRR"]])])))*MAX_Nsci$BIOM_Nsci
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Nsci$RESP<-0.01*2.3 * exp(0.0693*PARAMS[["SST"]][2:nrow(PARAMS[["SST"]]),"seasonal_AVG"]) * MIN_Nsci$BIOM_Nsci
MAX_Nsci$RESP-MIN_Nsci$RESP 
#Limits of egestion:
MIN_Nsci$EGEST<-0.1*MAX_Nsci$ING #checked for highest possible ingestion
MAX_Nsci$EGEST<-MAX_Nsci$RESP
MAX_Nsci$EGEST-MIN_Nsci$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Nsci$ING-MIN_Nsci$RESP-MIN_Nsci$EGEST
#NOTES:
#ACCEPTED
#but fraction of POC ingestion should be specified bc else eats up entire primary production
#--> >0.6*ING (arbitrary)

#####################################
####COPEPODS
#####################################
#-->Production
#-->Ingestion
#-->Specific ingestion of diatoms & microzooplankton
#-->Respiration
paramsCop<-c("BIOM_Cop","PROD","ING","ING_Mic","ING_Dia","RESP")
MIN_Cop<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsCop)+1),
                               dimnames=list(NULL,c("season",paramsCop))))
MIN_Cop$season<-seasons
MAX_Cop<-MIN_Cop
MAX_Cop[,"BIOM_Cop"]<-MIN_Cop[,"BIOM_Cop"]<-PARAMS[["CBIOM"]][,"Cop"]
MAX_Cop[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Cop"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Cop$BIOM_Cop
MAX_Cop[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Cop"),
                                                     2:ncol(PARAMS[["MIR"]])])))*MAX_Cop$BIOM_Cop
MAX_Cop[,"ING_Mic"]<-PARAMS[["SpecING_Cop"]][,"maxCR_Mic {m3/mgC copepod/day}"]*
                     PARAMS[["SpecING_Cop"]][,"PREY_C.Mic {mg C/m3}"]*MAX_Cop$BIOM_Cop
MIN_Cop[,"ING_Mic"]<-PARAMS[["SpecING_Cop"]][,"minCR_Mic {m3/mgC copepod/day}"]*
                     PARAMS[["SpecING_Cop"]][,"PREY_C.Mic {mg C/m3}"]*MIN_Cop$BIOM_Cop
MAX_Cop[,"ING_Dia"]<-PARAMS[["SpecING_Cop"]][,"maxCR_Dia {m3/mgC copepod/day}"]*
  PARAMS[["SpecING_Cop"]][,"PREY_C.Dia {mg C/m3}"]*MAX_Cop$BIOM_Cop
MIN_Cop[,"ING_Dia"]<-PARAMS[["SpecING_Cop"]][,"minCR_Dia {m3/mgC copepod/day}"]*
  PARAMS[["SpecING_Cop"]][,"PREY_C.Dia {mg C/m3}"]*MIN_Cop$BIOM_Cop
#Does minimum ingestion of Dia+Mic fit the maximum ingestion of copepods?
MAX_Cop$ING-MIN_Cop$ING_Mic-MIN_Cop$ING_Dia
#YES, except for spring 2010 (=huge diatom bloom)
#####SOLUTION#####
#Reduce diatom and microzooplankton concentration by 70% for spring 2010:
DiaSpring2010<-0.30 *
 PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),"PREY_C.Dia {mg C/m3}"]
MicSpring2010<-0.30 *
  PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),"PREY_C.Mic {mg C/m3}"]
MAX_Cop[which(MAX_Cop[,"season"]=="2010 spring"),"ING_Dia"]<-
  PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),
              "maxCR_Dia {m3/mgC copepod/day}"]*DiaSpring2010*
            MAX_Cop[which(MAX_Cop[,"season"]=="2010 spring"),"BIOM_Cop"]
MIN_Cop[which(MIN_Cop[,"season"]=="2010 spring"),"ING_Dia"]<-
  PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),
                    "minCR_Dia {m3/mgC copepod/day}"]*DiaSpring2010*
                  MIN_Cop[which(MIN_Cop[,"season"]=="2010 spring"),"BIOM_Cop"]
MAX_Cop[which(MAX_Cop[,"season"]=="2010 spring"),"ING_Mic"]<-
  PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),
                      "maxCR_Mic {m3/mgC copepod/day}"]*MicSpring2010*
                  MAX_Cop[which(MAX_Cop[,"season"]=="2010 spring"),"BIOM_Cop"]
MIN_Cop[which(MIN_Cop[,"season"]=="2010 spring"),"ING_Mic"]<-
  PARAMS[["SpecING_Cop"]][which(PARAMS[["SpecING_Cop"]][,"season"]=="2010 spring"),
                        "minCR_Mic {m3/mgC copepod/day}"]*MicSpring2010*
                MIN_Cop[which(MIN_Cop[,"season"]=="2010 spring"),"BIOM_Cop"]

MAX_Cop[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Cop"),
                                                      2:ncol(PARAMS[["MRR"]])])))*MAX_Cop$BIOM_Cop
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Cop$RESP<-0.01*2.3 * exp(0.0693*PARAMS[["SST"]][2:nrow(PARAMS[["SST"]]),"seasonal_AVG"]) * MIN_Cop$BIOM_Cop
MAX_Cop$RESP-MIN_Cop$RESP 
#Limits of egestion:
MIN_Cop$EGEST<-0.1*MAX_Cop$ING #checked for highest possible ingestion
MAX_Cop$EGEST<-MAX_Cop$RESP
MAX_Cop$EGEST-MIN_Cop$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Cop$ING-MIN_Cop$RESP-MIN_Cop$EGEST
#NOTES:
#ACCEPTED

#####################################
####TUNICATES
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsTun<-c("BIOM_Tun","PROD","ING","RESP")
MIN_Tun<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsTun)+1),
                               dimnames=list(NULL,c("season",paramsTun))))
MIN_Tun$season<-seasons
MAX_Tun<-MIN_Tun
MAX_Tun[,"BIOM_Tun"]<-MIN_Tun[,"BIOM_Tun"]<-PARAMS[["CBIOM"]][,"Tun"]
MAX_Tun[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Tun"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Tun$BIOM_Tun
MAX_Tun[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Tun"),
                                                     2:ncol(PARAMS[["MIR"]])])))*MAX_Tun$BIOM_Tun
MAX_Tun[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Tun"),
                                                      2:ncol(PARAMS[["MRR"]])])))*MAX_Tun$BIOM_Tun
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Tun$RESP<-0.01*2.3 * exp(0.0693*PARAMS[["SST"]][2:nrow(PARAMS[["SST"]]),"seasonal_AVG"]) * MIN_Tun$BIOM_Tun
MAX_Tun$RESP-MIN_Tun$RESP 
#Limits of egestion:
MIN_Tun$EGEST<-1/3*MAX_Tun$RESP #checked for highest possible respiration
MAX_Tun$EGEST<-MAX_Tun$RESP
MAX_Tun$EGEST-MIN_Tun$EGEST
#PROBLEM: min EGEST > max EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Tun$ING-MIN_Tun$RESP-MIN_Tun$EGEST
#NOTES:
#ACCEPTED

#####################################
####CLADOCERANS
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsClado<-c("BIOM_Clado","PROD","ING","RESP")
MIN_Clado<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsClado)+1),
                               dimnames=list(NULL,c("season",paramsClado))))
MIN_Clado$season<-seasons
MAX_Clado<-MIN_Clado
MAX_Clado[,"BIOM_Clado"]<-MIN_Clado[,"BIOM_Clado"]<-PARAMS[["CBIOM"]][,"Clado"]
MAX_Clado[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Clado"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Clado$BIOM_Clado
MAX_Clado[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Clado"),
                                                     2:ncol(PARAMS[["MIR"]])])))*MAX_Clado$BIOM_Clado
MAX_Clado[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Clado"),
                                                      2:ncol(PARAMS[["MRR"]])])))*MAX_Clado$BIOM_Clado
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Clado$RESP<-0.01*2.3 * exp(0.0693*PARAMS[["SST"]][2:nrow(PARAMS[["SST"]]),"seasonal_AVG"]) * MIN_Clado$BIOM_Clado
MAX_Clado$RESP-MIN_Clado$RESP 
#Limits of egestion:
MIN_Clado$EGEST<-0.1*MAX_Clado$ING #checked for highest possible ingestion
MAX_Clado$EGEST<-MAX_Clado$RESP
MAX_Clado$EGEST-MIN_Clado$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Clado$ING-MIN_Clado$RESP-MIN_Clado$EGEST
#NOTES:
#ACCEPTED

#####################################
####BIVALVE LARVAE
#####################################
#-->Production
#-->Ingestion
#-->Respiration
#(-->Specific Ingestion: not necessary to check as defined as fraction of ingestion)
paramsBiv<-c("BIOM_Biv","PROD","ING","SpecING_NanoEuk","SpecING_Phae","SpecING_Cil","SpecING_Bac","RESP")
MIN_Biv<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsBiv)+1),
                               dimnames=list(NULL,c("season",paramsBiv))))
MIN_Biv$season<-seasons
MAX_Biv<-MIN_Biv
MAX_Biv[,"BIOM_Biv"]<-MIN_Biv[,"BIOM_Biv"]<-PARAMS[["CBIOM"]][,"Biv"]
MAX_Biv[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Biv"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Biv$BIOM_Biv
MAX_Biv[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Biv"),
                                                     2:ncol(PARAMS[["MIR"]])])))*MAX_Biv$BIOM_Biv
MAX_Biv[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Biv"),
                                                     2:ncol(PARAMS[["MRR"]])])))*MAX_Biv$BIOM_Biv
#specific ingestions assuming maximum ingestion and minimum ingestion of nanoeukaryotes:
MIN_Biv[,"SpecING_NanoEuk"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Dia+Phae+Imp"),"MIN"]*MAX_Biv$ING
MIN_Biv[,"SpecING_Phae"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Phae"),"MIN"]*MIN_Biv$SpecING_NanoEuk
MAX_Biv[,"SpecING_Phae"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Phae"),"MAX"]*MIN_Biv$SpecING_NanoEuk
MIN_Biv[,"SpecING_Cil"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Cil"),"MIN"]*MAX_Biv$ING
MAX_Biv[,"SpecING_Cil"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Cil"),"MAX"]*MAX_Biv$ING
MAX_Biv[,"SpecING_Bac"]<-
  PARAMS[["SpecING_Biv"]][which(PARAMS[["SpecING_Biv"]][,"prey"]=="Bac"),"MAX"]*MAX_Biv$ING
#Do specific ingestion of ciliates fit with prey production available for grazing?
MAX_Cil$ING-MIN_Cil$RESP-MIN_Cil$EGEST-MIN_Biv$SpecING_Cil
#Lower bound of respiration-Vézina & Savenkoff (1999):
#<<At least 20% of ingestion>>
MIN_Biv$RESP<-0.2*MAX_Biv$ING
MAX_Biv$RESP-MIN_Biv$RESP
#-->DOES NOT FIT BUT IN THAT CASE BIVALVE INGESTION WILL NEVER REACH MIR
#Limits of egestion:
MIN_Biv$EGEST<-0.1*MAX_Biv$ING #checked for highest possible ingestion
MAX_Biv$EGEST<-MAX_Biv$RESP
MAX_Biv$EGEST-MIN_Biv$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Biv$ING-MIN_Biv$RESP-MIN_Biv$EGEST
#NOTES:
#ACCEPTED
#although ingestion too high for respiration therefore negative values: 
#Not a problem bc lower egestion/respiration limits are defined as function of ingestion
#=ingestion will never reach MIR

#####################################
####GASTROPOD LARVAE
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsGastr<-c("BIOM_Gastr","PROD","ING","RESP")
MIN_Gastr<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsGastr)+1),
                                dimnames=list(NULL,c("season",paramsGastr))))
MIN_Gastr$season<-seasons
MAX_Gastr<-MIN_Gastr
MAX_Gastr[,"BIOM_Gastr"]<-MIN_Gastr[,"BIOM_Gastr"]<-PARAMS[["CBIOM"]][,"Gastr"]
MAX_Gastr[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Gastr"),
                                                       2:ncol(PARAMS[["MGR"]])])))*MAX_Gastr$BIOM_Gastr
MAX_Gastr[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Gastr"),
                                                      2:ncol(PARAMS[["MIR"]])])))*MAX_Gastr$BIOM_Gastr
MAX_Gastr[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Gastr"),
                                                       2:ncol(PARAMS[["MRR"]])])))*MAX_Gastr$BIOM_Gastr
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Gastr$RESP<-0.2*MAX_Gastr$ING
MAX_Gastr$RESP-MIN_Gastr$RESP 
#Limits of egestion:
MIN_Gastr$EGEST<-0.1*MAX_Gastr$ING #checked for highest possible ingestion
MAX_Gastr$EGEST<-MAX_Gastr$RESP
MAX_Gastr$EGEST-MIN_Gastr$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Gastr$ING-MIN_Gastr$RESP-MIN_Gastr$EGEST
#NOTES:
#ACCEPTED

#####################################
####POLYCHAETE LARVAE
#####################################
#-->Production
#-->Ingestion
#-->Respiration
paramsPoly<-c("BIOM_Poly","PROD","ING","RESP")
MIN_Poly<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsPoly)+1),
                                dimnames=list(NULL,c("season",paramsPoly))))
MIN_Poly$season<-seasons
MAX_Poly<-MIN_Poly
MAX_Poly[,"BIOM_Poly"]<-MIN_Poly[,"BIOM_Poly"]<-PARAMS[["CBIOM"]][,"Poly"]
MAX_Poly[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Poly"),
                                                        2:ncol(PARAMS[["MGR"]])])))*MAX_Poly$BIOM_Poly
MAX_Poly[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Poly"),
                                                      2:ncol(PARAMS[["MIR"]])])))*MAX_Poly$BIOM_Poly
#MRR defined as average of Biv and Gastr MRR:
MRR_Poly<-data.frame(Biv=as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Biv"),
                                                            2:ncol(PARAMS[["MRR"]])]))),
                     Gastr=as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Gastr"),
                                                              2:ncol(PARAMS[["MRR"]])]))))
MAX_Poly[,"RESP"]<-apply(MRR_Poly,1,mean)*MAX_Poly$BIOM_Poly
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Poly$RESP<-0.2*MAX_Poly$ING
MAX_Poly$RESP-MIN_Poly$RESP 
#Limits of egestion:
MIN_Poly$EGEST<-0.1*MAX_Poly$ING #checked for highest possible ingestion
MAX_Poly$EGEST<-MAX_Poly$RESP
MAX_Poly$EGEST-MIN_Poly$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Poly$ING-MIN_Poly$RESP-MIN_Poly$EGEST
#NOTES:
#ACCEPTED
#although ingestion too high for respiration therefore negative values: 
#Not a problem bc lower egestion/respiration limits are defined as function of ingestion
#=ingestion will never reach MIR

#####################################
####HYDROMEDUSAE
#####################################
#-->Production
#-->Ingestion (MCR*prey concentration)
#-->Respiration
paramsHydro<-c("BIOM_Hydro","PROD","ING","RESP")
MIN_Hydro<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsHydro)+1),
                              dimnames=list(NULL,c("season",paramsHydro))))
MIN_Hydro$season<-seasons
MAX_Hydro<-MIN_Hydro
MAX_Hydro[,"BIOM_Hydro"]<-MIN_Hydro[,"BIOM_Hydro"]<-PARAMS[["CBIOM"]][,"Hydro"]
MAX_Hydro[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Hydro"),
                                                     2:ncol(PARAMS[["MGR"]])])))*MAX_Hydro$BIOM_Hydro
MAX_Hydro[,"ING"]<-as.vector(unname(t(PARAMS[["MCR"]][which(PARAMS[["MCR"]][,"Comp"]=="Hydro"),
                                            2:ncol(PARAMS[["MCR"]])])))*MAX_Hydro$BIOM_Hydro*
                   #as.vector(unname(t(PARAMS[["PREY_C"]][which(PARAMS[["PREY_C"]][,"Compartment"]=="Hydro"),
                                  #  2:(ncol(PARAMS[["PREY_C"]])-1)])))
                   as.vector(unname(t(PARAMS[["PREY_C"]][which(PARAMS[["PREY_C"]][,"Compartment"]=="Hydro"),
                                      "MEAN"])))
#Diatom biomass in spring 2010 is very high-->leads to huge ingestion compared to other seasons
#SOLUTION: use median diatom biomass in respective season
MAX_Hydro[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Hydro"),
                                                     2:ncol(PARAMS[["MRR"]])])))*MAX_Hydro$BIOM_Hydro
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Hydro$RESP<-0.2*MAX_Hydro$ING
MAX_Hydro$RESP-MIN_Hydro$RESP
#Limits of egestion:
MIN_Hydro$EGEST<-0.1*MAX_Hydro$ING #checked for highest possible ingestion
MAX_Hydro$EGEST<-MAX_Hydro$RESP
MAX_Hydro$EGEST-MIN_Hydro$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Hydro$ING-MIN_Hydro$RESP-MIN_Hydro$EGEST
#NOTES:
#ACCEPTED

#####################################
####M.LEIDYI
#####################################
#-->Production
#-->Ingestion (MCR*prey concentration)
#-->Respiration
#--> Diet
paramsMlei<-c("BIOM_Mlei","PROD","ING",paste("SpecING_",c("Dia","Tun","Nsci","Clado",
            "Cop","Biv","Gastr","Poly","Hydro"),sep=""),"RESP")
MIN_Mlei<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsMlei)+1),
                                dimnames=list(NULL,c("season",paramsMlei))))
MIN_Mlei$season<-seasons
MAX_Mlei<-MIN_Mlei
MAX_Mlei[,"BIOM_Mlei"]<-MIN_Mlei[,"BIOM_Mlei"]<-PARAMS[["CBIOM"]][,"Mlei"]
MAX_Mlei[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Mlei"),
                                                       2:ncol(PARAMS[["MGR"]])])))*MAX_Mlei$BIOM_Mlei
MAX_Mlei[,"ING"]<-as.vector(unname(t(PARAMS[["MCR"]][which(PARAMS[["MCR"]][,"Comp"]=="Mlei"),
                                                      2:ncol(PARAMS[["MCR"]])])))*MAX_Mlei$BIOM_Mlei*
                  as.vector(unname(t(PARAMS[["PREY_C"]][which(PARAMS[["PREY_C"]][,"Compartment"]=="Mlei"),
                                          "MEAN"])))
#Specific ingestion assuming maximum ingestion:
MAX_Mlei[which(MAX_Mlei[,"season"]=="2009 summer"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="MAX summer"),
                                     2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
                             MAX_Mlei[which(MAX_Mlei[,"season"]=="2009 summer"),"ING"]
MAX_Mlei[which(MAX_Mlei[,"season"]=="2010 summer"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="MAX summer"),
                                           2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
  MAX_Mlei[which(MAX_Mlei[,"season"]=="2010 summer"),"ING"]
MIN_Mlei[which(MIN_Mlei[,"season"]=="2009 summer"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="MIN summer"),
                                           2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
  MAX_Mlei[which(MIN_Mlei[,"season"]=="2009 summer"),"ING"]
MIN_Mlei[which(MIN_Mlei[,"season"]=="2010 summer"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="MIN summer"),
                                           2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
  MAX_Mlei[which(MIN_Mlei[,"season"]=="2010 summer"),"ING"]
MIN_Mlei[which(MIN_Mlei[,"season"]=="2009 autumn"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="autumn"),
                                           2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
  MAX_Mlei[which(MIN_Mlei[,"season"]=="2009 autumn"),"ING"]
MIN_Mlei[which(MIN_Mlei[,"season"]=="2010 autumn"),5:13]<-
  as.vector(unname(t(PARAMS[["DIET_Mlei"]][which(PARAMS[["DIET_Mlei"]][,"season"]=="autumn"),
                                           2:(ncol(PARAMS[["DIET_Mlei"]])-1)])))*
  MAX_Mlei[which(MIN_Mlei[,"season"]=="2010 autumn"),"ING"]
#Does minimum specific ingestion coincide with prey production (checked except for Dia)?
MAX_Nsci$ING-MIN_Nsci$RESP-MIN_Nsci$EGEST-MIN_Mlei$SpecING_Nsci 
MAX_Tun$ING-MIN_Tun$RESP-MIN_Tun$EGEST-MIN_Mlei$SpecING_Tun
MAX_Clado$ING-MIN_Clado$RESP-MIN_Clado$EGEST-MIN_Mlei$SpecING_Clado#absent in autumn 2010
MAX_Cop$ING-MIN_Cop$RESP-MIN_Cop$EGEST-MIN_Mlei$SpecING_Cop
MAX_Biv$ING-MIN_Biv$RESP-MIN_Biv$EGEST-MIN_Mlei$SpecING_Biv
MAX_Gastr$ING-MIN_Gastr$RESP-MIN_Gastr$EGEST-MIN_Mlei$SpecING_Gastr
MAX_Poly$ING-MIN_Poly$RESP-MIN_Poly$EGEST-MIN_Mlei$SpecING_Poly
MAX_Hydro$ING-MIN_Hydro$RESP-MIN_Hydro$EGEST-MIN_Mlei$SpecING_Hydro

MAX_Mlei[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Mlei"),
                                                       2:ncol(PARAMS[["MRR"]])])))*MAX_Mlei$BIOM_Mlei
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Mlei$RESP<-0.2*MAX_Mlei$ING
MAX_Mlei$RESP-MIN_Mlei$RESP
#-->fits except for spring 2010
#Limits of egestion:
MIN_Mlei$EGEST<-0.1*MAX_Mlei$ING #checked for highest possible ingestion
MAX_Mlei$EGEST<-MAX_Mlei$RESP
MAX_Mlei$EGEST-MIN_Mlei$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Mlei$ING-MIN_Mlei$RESP-MIN_Mlei$EGEST
#NOTES:
#ACCEPTED

#####################################
####P.PILEUS
#####################################
#-->Production
#-->Ingestion (MCR*prey concentration)
#-->Respiration
paramsPpil<-c("BIOM_Ppil","PROD","ING","RESP")
MIN_Ppil<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsPpil)+1),
                               dimnames=list(NULL,c("season",paramsPpil))))
MIN_Ppil$season<-seasons
MAX_Ppil<-MIN_Ppil
MAX_Ppil[,"BIOM_Ppil"]<-MIN_Ppil[,"BIOM_Ppil"]<-PARAMS[["CBIOM"]][,"Ppil"]
MAX_Ppil[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Ppil"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Ppil$BIOM_Ppil
MAX_Ppil[,"ING"]<-as.vector(unname(t(PARAMS[["MCR"]][which(PARAMS[["MCR"]][,"Comp"]=="Ppil"),
                                                     2:ncol(PARAMS[["MCR"]])])))*MAX_Ppil$BIOM_Ppil*
  as.vector(unname(t(PARAMS[["PREY_C"]][which(PARAMS[["PREY_C"]][,"Compartment"]=="Ppil"),
                                        "MEAN"])))
MAX_Ppil[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Ppil"),
                                                      2:ncol(PARAMS[["MRR"]])])))*MAX_Ppil$BIOM_Ppil
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Ppil$RESP<-0.2*MAX_Ppil$ING
MAX_Ppil$RESP-MIN_Ppil$RESP
#-->fits except for spring 2010
#Limits of egestion:
MIN_Ppil$EGEST<-0.1*MAX_Ppil$ING #checked for highest possible ingestion
MAX_Ppil$EGEST<-MAX_Ppil$RESP
MAX_Ppil$EGEST-MIN_Ppil$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Ppil$ING-MIN_Ppil$RESP-MIN_Ppil$EGEST
#NOTES:
#ACCEPTED

#####################################
####B.CUCUMIS
#####################################
#-->Production
#-->Ingestion (MCR*prey concentration)
#-->Respiration
paramsBcu<-c("BIOM_Bcu","PROD","ING","RESP")
MIN_Bcu<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsBcu)+1),
                               dimnames=list(NULL,c("season",paramsBcu))))
MIN_Bcu$season<-seasons
MAX_Bcu<-MIN_Bcu
MAX_Bcu[,"BIOM_Bcu"]<-MIN_Bcu[,"BIOM_Bcu"]<-PARAMS[["CBIOM"]][,"Bcu"]
MAX_Bcu[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Bcu"),
                                                      2:ncol(PARAMS[["MGR"]])])))*MAX_Bcu$BIOM_Bcu
MAX_Bcu[,"ING"]<-as.vector(unname(t(PARAMS[["MCR"]][which(PARAMS[["MCR"]][,"Comp"]=="Bcu"),
                                                     2:ncol(PARAMS[["MCR"]])])))*MAX_Bcu$BIOM_Bcu*
  as.vector(unname(t(PARAMS[["PREY_C"]][which(PARAMS[["PREY_C"]][,"Compartment"]=="Bcu"),
                                        "MEAN"])))
MAX_Bcu[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Bcu"),
                                                      2:ncol(PARAMS[["MRR"]])])))*MAX_Bcu$BIOM_Bcu
#Lower bound of respiration-Vézina & Platt (1988):
MIN_Bcu$RESP<-0.2*MAX_Bcu$ING
MAX_Bcu$RESP-MIN_Bcu$RESP
#-->fits except for spring 2010
#Limits of egestion:
MIN_Bcu$EGEST<-0.1*MAX_Bcu$ING #checked for highest possible ingestion
MAX_Bcu$EGEST<-MAX_Bcu$RESP
MAX_Bcu$EGEST-MIN_Bcu$EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
MAX_Bcu$ING-MIN_Bcu$RESP-MIN_Bcu$EGEST
#NOTES:
#ACCEPTED

#####################################
####HERRING
#####################################
#-->Production
#-->Ingestion
#-->Respiration
#-->Diet
paramsHer<-c("BIOM_Her","PROD","ING",paste("SpecING_",c("Dia","Tun","Clado",
                        "Cop","Biv","Gastr","Poly","Hydro"),sep=""),"RESP")
MIN_Her<-as.data.frame(matrix(NA,nrow=length(seasons),ncol=(length(paramsHer)+1),
                              dimnames=list(NULL,c("season",paramsHer))))
MIN_Her$season<-seasons
MAX_Her<-MIN_Her
MAX_Her[,"BIOM_Her"]<-MIN_Her[,"BIOM_Her"]<-PARAMS[["CBIOM"]][,"Her"]
MAX_Her[,"PROD"]<-as.vector(unname(t(PARAMS[["MGR"]][which(PARAMS[["MGR"]][,"Comp"]=="Her"),
                                                     2:ncol(PARAMS[["MGR"]])])))*MAX_Her$BIOM_Her
MAX_Her[,"ING"]<-as.vector(unname(t(PARAMS[["MIR"]][which(PARAMS[["MIR"]][,"Comp"]=="Her"),
                                                    2:ncol(PARAMS[["MIR"]])])))*MAX_Her$BIOM_Her
#Specific ingestion assuming maximum ingestion:
MAX_Her[which(MAX_Her[,"season"]=="2009 summer"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="MAX summer"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
                           MAX_Her[which(MAX_Her[,"season"]=="2009 summer"),"ING"]
MAX_Her[which(MAX_Her[,"season"]=="2010 summer"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="MAX summer"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
  MAX_Her[which(MAX_Her[,"season"]=="2010 summer"),"ING"]
MIN_Her[which(MIN_Her[,"season"]=="2009 summer"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="MIN summer"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
  MAX_Her[which(MIN_Her[,"season"]=="2009 summer"),"ING"]
MIN_Her[which(MIN_Mlei[,"season"]=="2010 summer"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="MIN summer"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
  MAX_Her[which(MIN_Her[,"season"]=="2010 summer"),"ING"]
MIN_Her[which(MIN_Her[,"season"]=="2009 autumn"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="autumn"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
  MAX_Her[which(MIN_Her[,"season"]=="2009 autumn"),"ING"]
MIN_Her[which(MIN_Her[,"season"]=="2010 autumn"),5:12]<-
  as.vector(unname(t(PARAMS[["DIET_Her"]][which(PARAMS[["DIET_Her"]][,"season"]=="autumn"),
                                           2:(ncol(PARAMS[["DIET_Her"]])-1)])))*
  MAX_Her[which(MIN_Her[,"season"]=="2010 autumn"),"ING"]
#Does minimum specific ingestion coincide with prey production (checked except for Dia)?
MAX_Tun$ING-MIN_Tun$RESP-MIN_Tun$EGEST-MIN_Her$SpecING_Tun
MAX_Clado$ING-MIN_Clado$RESP-MIN_Clado$EGEST-MIN_Her$SpecING_Clado#absent in autumn 2010
MAX_Cop$ING-MIN_Cop$RESP-MIN_Cop$EGEST-MIN_Her$SpecING_Cop
MAX_Biv$ING-MIN_Biv$RESP-MIN_Biv$EGEST-MIN_Her$SpecING_Biv
MAX_Gastr$ING-MIN_Gastr$RESP-MIN_Gastr$EGEST-MIN_Her$SpecING_Gastr
MAX_Poly$ING-MIN_Poly$RESP-MIN_Poly$EGEST-MIN_Her$SpecING_Poly
MAX_Hydro$ING-MIN_Hydro$RESP-MIN_Hydro$EGEST-MIN_Her$SpecING_Hydro


MAX_Her[,"RESP"]<-as.vector(unname(t(PARAMS[["MRR"]][which(PARAMS[["MRR"]][,"Comp"]=="Her"),
                                                     2:ncol(PARAMS[["MRR"]])])))*MAX_Her$BIOM_Her
#Lower bound of respiration-NO LITTERATURE SOURCE:

#Limits of egestion (Klumpp & von Westernhagen, 1986):
MIN_Her$EGEST<-0.066*MAX_Her$ING #checked for highest possible respiration
MAX_Her$EGEST<-MAX_Her$RESP
MAX_Her$EGEST-MIN_Her$EGEST
#PROBLEM: min EGEST > max EGEST
#What production is available for grazing assuming max. ingestion & min. respiration+egestion
#MAX_Tun$ING-MIN_Tun$RESP-MIN_Tun$EGEST
#NOTES:
#ACCEPTED
#but lack of lower boundary for respiration