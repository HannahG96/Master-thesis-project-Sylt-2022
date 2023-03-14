################# IN SITU RATES ##############################

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/")

#load packages
library(varhandle)
library(writexl)
library(data.table)

#########################################################
###DOC EXSUDATION RATE OF PHYTOPLANKTON {fraction of NPP}
######################################## Tillmann et al. (2000)
#Get carbon contribution of Phaeocystis sp. to total phytoplankton across spring 2009-winter 2010/11
source("Phyto-_Microplankton.spp_J.Rick/format_PhytoMicroplanktonData.R")
rm.all.but(keep=c("TOTCARB.SRB","CARBcontr.MICRO_SRB"), 
           envir=.GlobalEnv, keep_functions=FALSE, gc_limit=100,regex="auto")
seasons<-c("2009 spring","2009 summer","2009 autumn","2010 winter",
           "2010 spring","2010 summer","2010 autumn","2011 winter")
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
    PHAE<-TOTCARB.SRB[which(TOTCARB.SRB$season==seasons[i]),
            c("year","season","month","date","Phaeocystis_sp._C","Phytoplankton")]
  }else{
    TOTCARB<-rbind(TOTCARB,TOTCARB.SRB[which(TOTCARB.SRB$season==seasons[i]),])
    PHAE<-rbind(PHAE,TOTCARB.SRB[which(TOTCARB.SRB$season==seasons[i]),
                      c("year","season","month","date","Phaeocystis_sp._C","Phytoplankton")])}}

#Calculate C contribution of Phaeocystis to total phytoplankton at each sampling date:
PHAE$frac_PHAE<-PHAE$Phaeocystis_sp._C/PHAE$Phytoplankton

#Annotate the corresponding DOC exsudation rate at each sampling date
# 0.59-0.74 if Phae c/Phyto C >= 0.75
# 0.08-0.27 if Phae c/Phyto C < 0.75
PHAE$maxDOCexs<-PHAE$minDOCexs<-NA
for(i in 1:nrow(PHAE)){
  if(PHAE[i,"frac_PHAE"]>=0.75){
    PHAE[i,"maxDOCexs"]<-0.74
    PHAE[i,"minDOCexs"]<-0.59
  }else{
    PHAE[i,"maxDOCexs"]<-0.27
    PHAE[i,"minDOCexs"]<-0.08 }}

#Calculate average min and max exsudation rate of each season:
#-->average rates are based on monthly averages
DOCexs<-data.frame(season=seasons,frac_PHAE=NA,frac_PHAE.sd=NA,minDOCexs=NA,maxDOCexs=NA)
for(i in 1:length(seasons)){
  mydat<-PHAE[which(PHAE$season==seasons[i]),]
  DOCexs[i,"frac_PHAE"]<-mean(unname(tapply(mydat$frac_PHAE,INDEX=mydat$month,FUN=mean)))
  DOCexs[i,"frac_PHAE.sd"]<-sd(unname(tapply(mydat$frac_PHAE,INDEX=mydat$month,FUN=mean)))
  DOCexs[i,"minDOCexs"]<-mean(unname(tapply(mydat$minDOCexs,INDEX=mydat$month,FUN=mean)))
  DOCexs[i,"maxDOCexs"]<-mean(unname(tapply(mydat$maxDOCexs,INDEX=mydat$month,FUN=mean)))}

#Export DOC exsudation rates as excel file:
#write_xlsx(DOCexs, 
         #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/DOCexs.xlsx")

###########################################################################
#### GPP AND GPP-NPP {mg C/m3/day}
######################################################### Asmus & Asmus (2016)
source("R functions/get_estimates_MULTIparam.R")
#Import seasonal average, min and max GPP/NPP estimates:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/In situ rates/")
GPP_NPP<-fread(file = "GPP_NPP.csv", na.strings = "", dec = "," , data.table = FALSE)
#Reorder rows chronologically:
for(i in 1:length(seasons)){
  if(i==1){
    orderGPP_NPP<-as.data.frame(GPP_NPP[which(GPP_NPP$season==seasons[i]),])
  }else{
    orderGPP_NPP<-rbind(orderGPP_NPP,GPP_NPP[which(GPP_NPP$season==seasons[i]),])}}
#Min/Max monthly average carbon contribution of diatoms & Phaeocystis to total phytoplankton 
#biomass (in model) across season:
fracPHY<-merge(PHAE,TOTCARB.SRB[,c("date","Diatoms")],by="date",all=FALSE)
fracPHY$frac_DIA<-fracPHY$Diatoms/fracPHY$Phytoplankton
fracPHY.avg<-get_estimates_MULTIparam(dat=fracPHY,timeseries=2009:2011,
                                  stations=c("ST1"),colnames.params=c("frac_PHAE","frac_DIA"))[[1]]
fracPHY.avg<-fracPHY.avg[3:26,c("year","month","frac_PHAE_monthly_AVG","frac_DIA_monthly_AVG")]
colnames(fracPHY.avg)<-c("year","month","frac_PHAE","frac_DIA")

#Export GPP/NPP and phytoplankton fraction estimates as excel file:
write_xlsx(orderGPP_NPP, 
  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/GPP_NPP.xlsx")
write_xlsx(fracPHY.avg, 
           path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/fracPHY.xlsx")
##########################################################################################
###GRAZING OF PHYTOPLANKTON BY MICROZOOPLANKTON {fraction of phytoplankton standing stock}
################################################ Loebl & Van Beusekom (2008)    
#Explore the amount of grazed phytoplankton according to the min/max grazing rates defined
phyGRAZING<-TOTCARB[,c("season","month","date","Phytoplankton")]
phyGRAZING$maxGRAZING<-phyGRAZING$minGRAZING<-phyGRAZING$meanGRAZING<-NA
phyGRAZING[c(which(phyGRAZING$month==3),which(phyGRAZING$month==4),which(phyGRAZING$month==5)),
           "maxGRAZING"]<-0.66
phyGRAZING[c(which(phyGRAZING$month==3),which(phyGRAZING$month==4),which(phyGRAZING$month==5)),
           "minGRAZING"]<-0.01
phyGRAZING[c(which(phyGRAZING$month==3),which(phyGRAZING$month==4),which(phyGRAZING$month==5)),
           "meanGRAZING"]<-(0.01+0.04/3+0.66)/3
phyGRAZING[c(which(phyGRAZING$month==6),which(phyGRAZING$month==7),which(phyGRAZING$month==8)),
           "maxGRAZING"]<-1.23
phyGRAZING[c(which(phyGRAZING$month==6),which(phyGRAZING$month==7),which(phyGRAZING$month==8)),
           "minGRAZING"]<-0.51
phyGRAZING[c(which(phyGRAZING$month==6),which(phyGRAZING$month==7),which(phyGRAZING$month==8)),
           "meanGRAZING"]<-(1.23+0.51)/2
phyGRAZING[c(which(phyGRAZING$month==9),which(phyGRAZING$month==10),which(phyGRAZING$month==11)),
           "meanGRAZING"]<-0.17
phyGRAZING$maxPhyGRAZED<-phyGRAZING$Phytoplankton*phyGRAZING$maxGRAZING
phyGRAZING$minPhyGRAZED<-phyGRAZING$Phytoplankton*phyGRAZING$minGRAZING
phyGRAZING$meanPhyGRAZED<-phyGRAZING$Phytoplankton*phyGRAZING$meanGRAZING

#Calculate mean seasonal grazing rates (based on monthly averages):
PhyGRAZING<-data.frame(season=seasons,meanGRAZING=NA,minGRAZING=NA,maxGRAZING=NA,
                       meanPhyGRAZED=NA,minPhyGRAZED=NA,maxPhyGRAZED=NA)
for(i in 1:length(seasons)){
  mydat<-phyGRAZING[which(phyGRAZING$season==seasons[i]),]
  PhyGRAZING[i,"meanGRAZING"]<-mean(unname(tapply(mydat$meanGRAZING,INDEX=mydat$month,FUN=mean)))
  PhyGRAZING[i,"maxGRAZING"]<-mean(unname(tapply(mydat$maxGRAZING,INDEX=mydat$month,FUN=mean)))
  PhyGRAZING[i,"minGRAZING"]<-mean(unname(tapply(mydat$minGRAZING,INDEX=mydat$month,FUN=mean)))
  PhyGRAZING[i,"meanPhyGRAZED"]<-mean(unname(tapply(mydat$meanPhyGRAZED,INDEX=mydat$month,FUN=mean)))
  PhyGRAZING[i,"maxPhyGRAZED"]<-mean(unname(tapply(mydat$maxPhyGRAZED,INDEX=mydat$month,FUN=mean)))
  PhyGRAZING[i,"minPhyGRAZED"]<-mean(unname(tapply(mydat$minPhyGRAZED,INDEX=mydat$month,FUN=mean)))}

#Export Grazing rates as excel file:
#write_xlsx(PhyGRAZING, 
             #path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/PhyGRAZING.xlsx")

#############################################################
###PELAGIC PRODUCTION vs. PELAGIC RESPIRATION {dimensionless}
##############################################Loebl et al. (2007)
PtoR<-TOTCARB[,c("season","month","date","Diatoms","Phaeocystis_sp._C","Phytoplankton")]

#Define P/R for each season (based on monthly average PtoR):
#-->Diatom/Phaeocystis concentrations > 100 mg C/m3 during spring are defined as bloom
PtoR$maxPtoR<-PtoR$minPtoR<-NA
for(i in 1:nrow(PtoR)){
if(PtoR[i,"season"]=="2009 spring" || PtoR[i,"season"]=="2010 spring"){
    if(PtoR[i,"Diatoms"]>100 && PtoR[i,"Phaeocystis_sp._C"]>100){
      PtoR[i,"maxPtoR"]<-(5.3+4.6)/2
      PtoR[i,"minPtoR"]<-(3.9+3)/2
    }else if(PtoR[i,"Diatoms"]>100){
      PtoR[i,"maxPtoR"]<-5.3
      PtoR[i,"minPtoR"]<-3.9
    }else if(PtoR[i,"Phaeocystis_sp._C"]>100){
      PtoR[i,"maxPtoR"]<-4.6
      PtoR[i,"minPtoR"]<-3.0
    }else{
      PtoR[i,"maxPtoR"]<-5.1
      PtoR[i,"minPtoR"]<-4.1}}
if(PtoR[i,"season"]=="2009 summer" || PtoR[i,"season"]=="2010 summer"){
  PtoR[i,"maxPtoR"]<-2.8
  PtoR[i,"minPtoR"]<-1.8}
if(PtoR[i,"season"]=="2009 autumn" || PtoR[i,"season"]=="2010 autumn"){
  PtoR[i,"maxPtoR"]<-1.9
  PtoR[i,"minPtoR"]<-0.6}
if(PtoR[i,"season"]=="2010 winter" || PtoR[i,"season"]=="2011 winter"){
  PtoR[i,"maxPtoR"]<-0.7
  PtoR[i,"minPtoR"]<-0.3}
}#END OF FOR LOOP
PToR<-data.frame(season=seasons,minPtoR=NA,maxPtoR=NA)
for(i in 1:length(seasons)){
  mydat<-PtoR[which(PtoR$season==seasons[i]),]
  PToR[i,"minPtoR"]<-mean(unname(tapply(mydat$minPtoR,INDEX=mydat$month,FUN=mean)))
  PToR[i,"maxPtoR"]<-mean(unname(tapply(mydat$maxPtoR,INDEX=mydat$month,FUN=mean)))}

#Export P/R-ratios as excel file:
#write_xlsx(PToR, 
 #path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/PToR.xlsx")

#####################################################################################
###BACTERIAL CARBON DEMAND (cf. DOC UPTAKE) {fraction of exsudated phytoplankton DOC}
############################################ Sintes et al. (2010)
BCD<-data.frame(year=c(rep(2009,10),rep(2010,12),rep(2011,2)),
                yearS=c(rep(2009,9),rep(2010,12),rep(2011,3)),month=c(3:12,1:12,1,2),
                season=rep(c(rep("spring",3),rep("summer",3),rep("autumn",3),rep("winter",3)),2),
                BCD=NA)
BCD$season<-paste(BCD$yearS,BCD$season,sep=" ")
BCD[which(BCD$month<=7),"BCD"]<-1/0.98
BCD[which(BCD$month>=8),"BCD"]<-1/0.48
DOCuptake<-data.frame(season=seasons,DOCuptake=NA,SD=NA)
for(i in 1:length(seasons)){
  DOCuptake[i,"DOCuptake"]<-mean(BCD[which(BCD[,"season"]==seasons[i]),"BCD"])
  DOCuptake[i,"SD"]<-sd(BCD[which(BCD[,"season"]==seasons[i]),"BCD"])}
#Export DOC uptake rates as excel file:
#write_xlsx(DOCuptake, 
  #path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/DOCuptake.xlsx")

#######################################################################
###BACTERIAL GROWTH EFFICIENCY {dimensionless}
#(=Bacterial Production/(Bacterial Production + Bacterial Respiration))
############################################ Sintes et al. (2010)
BacGrowthEff<-data.frame(season=seasons,minBGE=NA,maxBGE=NA)
for(i in 1:nrow(BacGrowthEff)){
  if(BacGrowthEff[i,"season"]=="2009 spring" || BacGrowthEff[i,"season"]=="2010 spring"){
    BacGrowthEff[i,"minBGE"]<- 0.16
    BacGrowthEff[i,"maxBGE"]<- 0.43
} else if(BacGrowthEff[i,"season"]=="2009 summer" || BacGrowthEff[i,"season"]=="2010 summer"){
    BacGrowthEff[i,"minBGE"]<- 0.1
    BacGrowthEff[i,"maxBGE"]<- 0.53
} else if(BacGrowthEff[i,"season"]=="2009 autumn" || BacGrowthEff[i,"season"]=="2010 autumn"){
    BacGrowthEff[i,"minBGE"]<- 0.06
    BacGrowthEff[i,"maxBGE"]<- 0.49
} else if(BacGrowthEff[i,"season"]=="2010 winter" || BacGrowthEff[i,"season"]=="2011 winter"){
    BacGrowthEff[i,"minBGE"]<- 0.06
    BacGrowthEff[i,"maxBGE"]<- 0.46}}

#Export Bacterial Growth efficiencies as excel file:
#write_xlsx(BacGrowthEff, 
      #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/BacGrowthEff.xlsx")

#######################################################################
###GROWTH RATE OF PHAEOCYSTIS SP. {per day}
############################################ Weisse & Scheffel-Möser (1990)
#-->growth rate is seasonally adapted based on an assumed Q10 of 1.93 (Tillmann, 2000)
#Import data of seasonal average SST:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/In situ rates/")
SST<-fread(file = "SST.csv", na.strings = "", dec = "," , data.table = FALSE)
SST<-SST[-1,]
#Calculate temperature-adapted min and max growth rates of Phaeocystis:
# Rt=Rt0*Q10^((t-t0)/10)
minGR<-0.033*24 #{per day}
maxGR<-0.098*24 #{per day}
minGR_season<-minGR*1.93^((SST[,"seasonal_AVG"]-9)/10)
maxGR_season<-maxGR*1.93^((SST[,"seasonal_AVG"]-9)/10)
GROWTH_Phae<-data.frame(season=seasons,MIN=minGR_season,MAX=maxGR_season)
#Export growth rates of Phaeocystis as Excel file:
#write_xlsx(GROWTH_Phae, 
 # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/GROWTH_Phae.xlsx")

#######################################################################
###GROWTH RATE OF DIATOMS {per day}
############################################ Stelfox-Widdicombe et al. (2004)
#-->growth rate is seasonally adapted based on an assumed Q10 of 1.93 (Tillmann, 2000)

#Calculate temperature-adapted min and max growth rates of diatoms:
# Rt=Rt0*Q10^((t-t0)/10)
minGR<-0.9 #{per day}
maxGR<-1.5 #{per day}
T0<-(7.1+8.9)/2
minGR_season<-minGR*1.93^((SST[,"seasonal_AVG"]-T0)/10)
maxGR_season<-maxGR*1.93^((SST[,"seasonal_AVG"]-T0)/10)
GROWTH_Dia<-data.frame(season=seasons,MIN=minGR_season,MAX=maxGR_season)
#Export growth rates of Phaeocystis as Excel file:
#write_xlsx(GROWTH_Dia, 
   # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/GROWTH_Dia.xlsx")

#######################################################################
###POC EXPORT {mg POC/m3/day}
############################################ Fonova et al. (2019); Postma (1981)
#<<The tide entering the Wadden Sea leaves around ~3.5*10^6 tons/year of suspended matter>>
#<<Estimates of the mean water volume entering the basin during flood and leaving during ebb
#(the tidal prism) vary from 4*10^8 to 6.3*10^8 m3>>

#Import POC data of SRB (indicating the deviation of the seasonal POC stock from the 
#annual POC stock):
POC<-fread(file = "POC.csv", na.strings = "", dec = "," , data.table = FALSE)

#Calculate the mean water volume leaving the bight per day:
dailyOUTFLOW<-(4+6.3)/2 * 10^8 * 2 #{m3/day}

#Calculate the mean amount of SPM transported out of the bight per day:
dailySPMEXPORT<-(3.5*10^6)/365 * 10^9 #{mg/day}

#Calculate the mean concentration of SPM transported out of the bight per day:
dailySPMEXPORTperm3<-dailySPMEXPORT/dailyOUTFLOW

#Convert daily SPM export into POC export based on linear regression of Spiekeroog 
#Backbarrier Tidal Flat data - Dellwig et al. (2006):
#POC {mg/L}~SPM{mg/L}
dailyPOCEXPORTperm3<-(0.1627051+0.0342496*(dailySPMEXPORTperm3*0.001))/0.001 #{mg/m3/day}

#Adjust daily POC export to seasons based on the deviation of seasonal POC stocks from
#annual POC stock:
POC$EXPORT<-POC$Dev_annualPOC*dailyPOCEXPORTperm3

#Export seasonal POC exports to the North Sea as Excel file:
#write_xlsx(POC, 
    #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/POC_EXPORTtoNorthSea.xlsx")
           