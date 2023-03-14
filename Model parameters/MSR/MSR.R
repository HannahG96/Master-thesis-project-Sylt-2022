##################################################################
####MAXIMUM SPECIFIC RATES of RR/IR/CR for model compartments#####
##################################################################

#!!Note: Typing error discovered, "R.octopuntata" vs."R.octopunctata" (??????)

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/MSR")

#load packages
library(data.table)
library(doBy)
library(ggplot2)
library(writexl)
library(varhandle)

#############################
#A)SPP-SPECIFIC CARBON WEIGHT

#Estimate avg carbon biomass of each species/func. group composing the compartments
CARB<-fread(file = "biomassCarbon.csv", na.strings = "", dec = "," , data.table = FALSE)#import df with length-carbon relationships
CARB$LWR<-sub(",",".",CARB$LWR,fixed=TRUE)
Unit<-apply(as.data.frame(CARB[,"[W]"]),1,FUN=function(x)strsplit(x,split="[",fixed=TRUE)[[1]][2])
Unit<-sub("]","",Unit,fixed=TRUE)
CBiom<-data.frame(Comp=CARB[,"Compartment"],Spp=CARB[,"Studied taxon"],Col=CARB[,"Column"],
                  biomassC=NA, Unit=Unit)

#Estimated species-specific carbon biomass based on LWR:
assumedL.hydro<-mean(c(9,4.5,5.5,6.405))#unknown hydromedusae size is estimated based on mean of known hydromedusae sizes
assumedDW.copepods<-mean(c(0.465,0.497,0.433,0.453))#for all copepod species for which no average C/W-ratio is available, it is assumed an average C/W-ratio of the compartment group

for(i in 1:nrow(CARB)){
  myspp<-CARB[i,"Column"]
  LWR.myspp<-parse(text=CARB[i,"LWR"])
  L<-CARB[i,"meanL"]
  if(is.na(L)==TRUE)L<-assumedL.hydro
  C.myspp<-eval(LWR.myspp)
  print(C.myspp)
  C.ratio<-CARB[i,"C/W-ratio"] #(if necessary) adjust biomass estimate based on DW/WW and C/WW ratios
  if(is.na(C.ratio)==TRUE && CARB[i,"Compartment"]=="Cop"){C.ratio<-0.15*assumedDW.copepods #DW/WW fraction of copepods is assumed 15%
  }else if(is.na(C.ratio)==TRUE){C.ratio<-1
  }else if(CARB[i,"Compartment"]=="Cop"){C.ratio<-0.15*C.ratio}
  print(C.ratio)
  CBiom[i,"biomassC"]<-C.myspp*C.ratio}
#Convert all carbon estimates in mg C
CBiom[which(CBiom[,"Unit"]=="g"),"biomassC"]<-CBiom[which(CBiom[,"Unit"]=="g"),"biomassC"]*1000#convert gram C in mg C
CBiom[which(CBiom[,"Unit"]=="ug"),"biomassC"]<-CBiom[which(CBiom[,"Unit"]=="ug"),"biomassC"]*0.001#convert ug C in mg C
CBiom[which(CBiom[,"Unit"]=="fg"),"biomassC"]<-CBiom[which(CBiom[,"Unit"]=="fg"),"biomassC"]*10^(-12)#convert fg C in mg C
CBiom<-CBiom[,1:4]
#Get average carbon biomass of Spionidae:
CBiom[which(CBiom[,"Col"]=="Spionidae"),"biomassC"]<-mean(CBiom[which(CBiom[,"Col"]=="Spionidae"),"biomassC"])
#Add assumed biomass of N. scintillans+ciliates
CBiom<-rbind(CBiom,c("Nsci","Noctiluca scintillans","N.scintillans","0.0002"))
CBiom<-rbind(CBiom,c("Cil","ciliates","ciliates","9.5e-06"))
CBiom$biomassC<-as.numeric(CBiom$biomassC)

#Export Specific carbon weight of species as Excelfile:
#write_xlsx(CBiom, 
       #    path="C:/Hannah/Biological Oceanography/Master Thesis/Project/METHODOLOGY/BodyCarb.xlsx")

################################################################
#B)COMPUTATION OF MSR FOR HET. COMPARTMENTS (except bacteria)

#Import df with MSR relationships
MSR<-fread(file = "MSR.csv", na.strings = "", dec = "," , data.table = FALSE)
colnames(MSR)[9]<-"MSR"
MSR[,"MSR"]<-sub(",",".",MSR[,"MSR"],fixed=TRUE)
MSR[,"Temperature"]<-sub(",",".",MSR[,"Temperature"],fixed=TRUE)
MSR$Temperature<-as.numeric(MSR$Temperature)
MSR[,"Conversion [mg C/mg C/day] or [m3/mg C/day]"]<-sub(",",".",MSR[,"Conversion [mg C/mg C/day] or [m3/mg C/day]"],fixed=TRUE)
MSR[,"Conversion [mg C/mg C/day] or [m3/mg C/day]"]<-sub("MSR","msr",MSR[,"Conversion [mg C/mg C/day] or [m3/mg C/day]"],fixed=TRUE)

#Create dfs to store species-specific MSR values
Umax<-as.data.frame(matrix(NA,nrow=0,ncol=6,
                        dimnames=list(NULL,c("Comp","Spp","Col","Weight","Umax","Temp"))))
Imax<-as.data.frame(matrix(NA,nrow=0,ncol=6,
                        dimnames=list(NULL,c("Comp","Spp","Col","Weight","Imax","Temp"))))
Cmax<-as.data.frame(matrix(NA,nrow=0,ncol=7,
                        dimnames=list(NULL,c("Comp","Spp","Col","Weight","Cmax","Temp","C_PREY"))))
Rmax<-as.data.frame(matrix(NA,nrow=0,ncol=6,
                        dimnames=list(NULL,c("Comp","Spp","Col","Weight","Rmax","Temp"))))

#Calculate weight-specific MSR of species:
comps<-c("Cil","Nsci","Cop","Clado","Tun","Poly","Biv","Gastr",
         "Hydro","Mlei","Ppil","Bcu","Her")
n<-1
for(i in 1:length(comps)){
  CBiom.mycomp<-as.data.frame(CBiom[which(CBiom[,"Comp"]==comps[i]),])

for(u in 1:nrow(CBiom.mycomp)){
  myspp<-CBiom.mycomp[u,"Col"]
  M<-as.numeric(CBiom.mycomp[u,"biomassC"])
  MSR.myspp<-as.data.frame(MSR[which(MSR[,"Column"]==myspp),])
  infos<-unname(as.vector(CBiom.mycomp[u,c("Comp","Spp","Col","biomassC")]))
  Umax[n,c("Comp","Spp","Col","Weight")]<-Imax[n,c("Comp","Spp","Col","Weight")]<-
  Cmax[n,c("Comp","Spp","Col","Weight")]<-Rmax[n,c("Comp","Spp","Col","Weight")]<-infos

#weight-specific growth rate:  
  msr<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="GR"),"MSR"]))#weight-specific MSR relationship
  conv.GR<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="GR"),"Conversion [mg C/mg C/day] or [m3/mg C/day]"]))#conversion to mg C/mg C/day
  if(is.na(conv.GR)==FALSE){Umax[n,"Umax"]<-conv.GR}else{Umax[n,"Umax"]<-msr}
 Umax[n,"Temp"]<-MSR.myspp[which(MSR.myspp[,"Rate R"]=="GR"),"Temperature"]

#weight-specific ingestion rate:
if(length(which(MSR.myspp[,"Rate R"]=="IR"))!=0){
  msr<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="IR"),"MSR"]))#weight-specific MSR relationship
  conv.IR<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="IR"),"Conversion [mg C/mg C/day] or [m3/mg C/day]"]))#conversion to mg C/mg C/day
  if(is.na(conv.IR)==FALSE){Imax[n,"Imax"]<-conv.IR}else{Imax[n,"Imax"]<-msr}
  Imax[n,"Temp"]<-MSR.myspp[which(MSR.myspp[,"Rate R"]=="IR"),"Temperature"]}

#weight-specific clearance rate:  
if(length(which(MSR.myspp[,"Rate R"]=="CR"))!=0){
  msr<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="CR"),"MSR"]))#weight-specific MSR relationship
  C_PREY<-1 #PREY CONCENTRATION IS SET 1 TO OBTAIN CLEARANCE RATE IN M3
  conv.CR<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="CR"),"Conversion [mg C/mg C/day] or [m3/mg C/day]"]))#conversion to mg C/mg C/day
  if(is.na(conv.CR)==FALSE){Cmax[n,"Cmax"]<-conv.CR}else{Cmax[n,"Cmax"]<-msr}
  Cmax[n,"C_PREY"]<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="CR"),"C_PREY [mg C/m3]"]))
  Cmax[n,"Temp"]<-MSR.myspp[which(MSR.myspp[,"Rate R"]=="CR"),"Temperature"]}

#weight-specific respiration rate:  
 msr<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="RR"),"MSR"]))#weight-specific MSR relationship
 conv.RR<-eval(parse(text=MSR.myspp[which(MSR.myspp[,"Rate R"]=="RR"),"Conversion [mg C/mg C/day] or [m3/mg C/day]"]))#conversion to mg C/mg C/day
 if(is.na(conv.RR)==FALSE){Rmax[n,"Rmax"]<-conv.RR}else{Rmax[n,"Rmax"]<-msr}
 Rmax[n,"Temp"]<-MSR.myspp[which(MSR.myspp[,"Rate R"]=="RR"),"Temperature"]
 n<-n+1
 print(n)
}}

#################################################################
###ADJUST MSR TO SEASONAL TEMPERATURE BASED ON Q10: spring 2009-winter 2010/11
#################################################################

#Import df with Q10 infos:
Q10<-fread(file = "Q10.csv", na.strings = "", dec = "," , data.table = FALSE)
Q10$Q10<-sub(",",".",Q10$Q10, fixed=TRUE)
Q10$Q10<-as.numeric(Q10$Q10)
#Import df with SST info:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/")
dat<-fread(file = "Langzeitdaten_J.Rick/Hydrochemistry/SRB_Hydrochemistry_Timeseries.csv", na.strings = "", dec = "," , data.table = FALSE)
#Make date column understandable for R
dat$Date<- as.Date (dat$Date , format = "%d.%m.%Y")
#Add year & month column
dat$year<-dat$month<-NA
for(i in 1:nrow(dat)){
  dat[i,"year"]<-strsplit(as.character(dat[i,"Date"]),split="-")[[1]][1]
  dat[i,"month"]<-strsplit(as.character(dat[i,"Date"]),split="-")[[1]][2]}
dat$year<-as.numeric(dat$year)
dat$month<-as.numeric(dat$month)
#Calculate seasonal SST estimates:
source("R functions/get_estimates_ONEparam.R")
data.SST<-get_estimates_ONEparam(dat=dat,timeseries=2009:2011,
                                 stations=c("Ferry Terminal","1","4"),colname.param="SST")[[2]]
data.SST$season<-paste(data.SST$year,data.SST$season,sep=" ")
seasonal_SST<-data.SST[c(2:9),c("year","season","seasonal_AVG","seasonal_SD","seasonal_pctSD")]
colnames(seasonal_SST)<-c("year","season","SST","SD","pctSD")

##############################
#Calculate seasonal MSRs

#1.Create function to temperature-adjust MSR based on Q10:
MSR_temp<-function(msr,q10,t0,t){
  msr.temp<-msr*q10^((t-t0)/10)
  return(msr.temp)}

#2.Create dfs to store temperature-adjusted MSR values:
seasons<-as.data.frame(matrix(NA,nrow=nrow(Umax),ncol=nrow(seasonal_SST),
                              dimnames=list(NULL,seasonal_SST[,"season"])))
Umax2009_10<-cbind(Umax,seasons)
Imax2009_10<-cbind(Imax,seasons)
Cmax2009_10<-cbind(Cmax,seasons)
Rmax2009_10<-cbind(Rmax,seasons)

#3. Calculate seasonal estimates for each species or compartment
for(i in 1:nrow(Umax)){
  mycomp<-Umax[i,"Comp"]
  myspp<-Umax[i,"Col"]
  q10.myspp<-Q10[which(Q10[,"Compartment"]==mycomp),]
  for(u in 1:nrow(seasonal_SST)){
    temp<-seasonal_SST[u,"SST"]
  #growth/production rate:
    myq10<-q10.myspp[which(q10.myspp[,"Rate"]=="PR"),"Q10"]
    t0<-Umax[i,"Temp"]
    R0<-Umax[i,"Umax"]
    Umax2009_10[i,(ncol(Umax)+u)]<-MSR_temp(msr=R0,q10=myq10,t0=t0,t=temp)
  #ingestion rate:
    myq10<-q10.myspp[which(q10.myspp[,"Rate"]=="IR"),"Q10"]
  if(length(myq10)!=0){  
    t0<-Imax[i,"Temp"]
    R0<-Imax[i,"Imax"]
    Imax2009_10[i,(ncol(Imax)+u)]<-MSR_temp(msr=R0,q10=myq10,t0=t0,t=temp)}
  #clearance rate:
    myq10<-q10.myspp[which(q10.myspp[,"Rate"]=="CR"),"Q10"]
  if(length(myq10)!=0){
    t0<-Cmax[i,"Temp"]
    R0<-Cmax[i,"Cmax"]
    Cmax2009_10[i,(ncol(Cmax)+u)]<-MSR_temp(msr=R0,q10=myq10,t0=t0,t=temp)}
  #respiration rate:
    myq10<-q10.myspp[which(q10.myspp[,"Rate"]=="RR"),"Q10"]
  if(length(myq10)!=0){
    t0<-Rmax[i,"Temp"]
    R0<-Rmax[i,"Rmax"]
    Rmax2009_10[i,(ncol(Rmax)+u)]<-MSR_temp(msr=R0,q10=myq10,t0=t0,t=temp)}}}

###################################################################
###CALCULATE BIOMASS-WEIGHTED SEASONAL MSR OF HET. COMPARTMENTS (except bacteria)
#-->compartment-msr as a weighted average of species-specific msr 
#   based on average biomass-contributions in each season
###################################################################

#1.Import infos about seasonal biomass contribution of mero-/zoo- 
#& gelatinous zooplankton:

source("Zooplankton_P.Martens/format_ZooplanktonData.R")#mesozooplankton
source("Langzeitdaten_J.Rick/Meroplankton/format_MeroplanktonData.R")#meroplankton
source("Langzeitdaten_J.Rick/Gelatinous Zooplankton/format_GelatinousZooData.R")#gelatinous zooplankton
rm.all.but(keep=c("seasonal_SST","Umax2009_10","Imax2009_10","Cmax2009_10","Rmax2009_10",
                  "CARBcontr.MERO","CARBcontr.ZOO","CARBcontr.GelZOO","comps"), 
           envir=.GlobalEnv, keep_functions=FALSE, gc_limit=100,regex="auto")
CARBcontr<-rbind(CARBcontr.MERO,CARBcontr.ZOO,CARBcontr.GelZOO)
for(i in 1:nrow(CARBcontr)){
  if(length(strsplit(CARBcontr[i,"Species"],split=" ",fixed=TRUE)[[1]])>1){
    print(CARBcontr[i,"Species"])
    rename<-paste(strsplit(CARBcontr[i,"Species"],split=" ",fixed=TRUE)[[1]][1],
                strsplit(CARBcontr[i,"Species"],split=" ",fixed=TRUE)[[1]][2],sep="_")
  rename<-sub(".","",rename,fixed=TRUE)
  CARBcontr[i,"Species"]<-rename
  print(CARBcontr[i,"Species"])}}

#comp.names<-data.frame(Comp=Umax2009_10$Comp,Species=Umax2009_10$Col)
#CARBcontr<-merge(comp.names[3:30,],CARBcontr,all=TRUE,by="Species")

#2.Create dfs to store biomass-weighted seasonal MSR of compartments:
UMAX<-IMAX<-CMAX<-RMAX<-
  cbind(data.frame(Comp=comps),
  as.data.frame(matrix(NA,nrow=length(comps), ncol=nrow(seasonal_SST),
        dimnames=list(NULL,seasonal_SST$season))))
#3.Calculate compartment-level GR/IR/CR/RR
mixed.comps<-c(rep(NA,2),"copepods","cladocerans",NA,"polychaetes",
               rep(NA,2),"Hydromedusae",rep(NA,4))
seasons<-seasonal_SST$season
for(i in 1:length(comps)){

  if(is.na(mixed.comps[i])==TRUE){
    UMAX[i,2:ncol(UMAX)]<-Umax2009_10[which(Umax2009_10[,"Comp"]==comps[i]),seasons]
    IMAX[i,2:ncol(IMAX)]<-Imax2009_10[which(Imax2009_10[,"Comp"]==comps[i]),seasons]
    CMAX[i,2:ncol(CMAX)]<-Cmax2009_10[which(Cmax2009_10[,"Comp"]==comps[i]),seasons]
    RMAX[i,2:ncol(RMAX)]<-Rmax2009_10[which(Rmax2009_10[,"Comp"]==comps[i]),seasons]
  }else{
fracs<-CARBcontr[which(CARBcontr[,"Compartment"]==mixed.comps[i]),]
umax<-Umax2009_10[which(Umax2009_10[,"Comp"]==comps[i]),]
imax<-Imax2009_10[which(Imax2009_10[,"Comp"]==comps[i]),]
cmax<-Cmax2009_10[which(Cmax2009_10[,"Comp"]==comps[i]),]
rmax<-Rmax2009_10[which(Rmax2009_10[,"Comp"]==comps[i]),]

 for(u in 1:nrow(seasonal_SST)){
fracs.season<-fracs[which(fracs[,"Season"]==seasons[u]),c("pct_Carbon","Species")]
fracs.season$Frac<-fracs.season$pct_Carbon/100
colnames(fracs.season)[2:3]<-c("Col","Frac")
#growth rate:
my.umax<-merge(umax[,c("Col",seasons[u])],fracs.season[,c("Col","Frac")],by="Col",all=TRUE)
#print(umax$Col)
umax.weighted<-my.umax[,seasons[u]]*my.umax[,"Frac"]
UMAX[i,seasons[u]]<-sum(umax.weighted) 
#ingestion rate:
my.imax<-merge(imax[,c("Col",seasons[u])],fracs.season[,c("Col","Frac")],by="Col")
print(imax$Col)
imax.weighted<-my.imax[,seasons[u]]*my.imax[,"Frac"]
IMAX[i,seasons[u]]<-sum(imax.weighted)
#clearance rate:
my.cmax<-merge(cmax[,c("Col",seasons[u])],fracs.season[,c("Col","Frac")],by="Col")
print(cmax$Col)
cmax.weighted<-my.cmax[,seasons[u]]*my.cmax[,"Frac"]
CMAX[i,seasons[u]]<-sum(cmax.weighted)
#respiration rate:
my.rmax<-merge(rmax[,c("Col",seasons[u])],fracs.season[,c("Col","Frac")],by="Col")
print(rmax$Col)
rmax.weighted<-my.rmax[,seasons[u]]*my.umax[,"Frac"]
RMAX[i,seasons[u]]<-sum(rmax.weighted)}}#END OF ELSE CLAUSE
}

#EXPORT SEASONAL MSR DATA FRAMES AS EXCEL FILES:
#UMAX:
#write_xlsx(UMAX, 
          # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/MGR.xlsx")
#IMAX:
#write_xlsx(IMAX, 
          # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/MIR.xlsx")
#CMAX:
#write_xlsx(CMAX, 
          # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/MCR.xlsx")
#UMAX:
#write_xlsx(RMAX, 
          # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/MRR.xlsx")


#####################################################################################
#CONTROLS FOR LOWER BOUNDARIES:
temp<-seasonal_SST$SST
#protozooplankton respiration
as.numeric(unname(RMAX[1,2:9]))-(0.01*6.4*exp(0.0693*temp))
#protozooplankton production
as.numeric(unname(UMAX[1,2:9]))-0.5*as.numeric(unname(IMAX[1,2:9]))
#meso-, meroplankton respiration
for(i in 2:8){a<-as.numeric(unname(RMAX[i,2:9]))-(0.01*2.3*exp(0.0693*temp))
print(RMAX[i,1])
print(a)}
#-->negative thresholds for meroplankton=define another lower boundary
#macrozooplankton respiration
#for(i in 9:12){a<-as.numeric(unname(RMAX[i,2:9]))-0.2*as.numeric(unname(CMAX[i,2:9]))
#print(RMAX[i,1])
#print(a)}
#macrozooplankton production
#for(i in 9:12){a<-as.numeric(unname(RMAX[i,2:9]))-as.numeric(0.2*unname(CMAX[i,2:9]))
#print(RMAX[i,1])
#print(a)}
as.numeric(CMAX[9,2:ncol(CMAX)])*12.80317
