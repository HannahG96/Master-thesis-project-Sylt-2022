###########################################################
#ANALYSE & FORMAT PELAGIC FISH TIME SERIES OF SRB
###########################################################
#-->sampled at List Reede (2007-2019): 23.01.2007-20.11.2019

#load packages
library(data.table)
library(ggplot2)
library(doBy)

#set working directory
working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/"
plot.working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Community composition/"
setwd(working_dir)

##############################################################################
###IMPORT & FORMAT DATA

##Herring (2007-2019)######################################################################
dat1<-fread(file = "Fish/data_herring.csv", na.strings = "", dec = "," , data.table = FALSE)
LWR<-fread(file = "Fish/LWR_herring.csv", na.strings = "", dec = "," , data.table = FALSE)#length-weight relationships for herring

dat1$date<-as.Date(apply(as.data.frame(dat1[,8]),1,FUN=function(x)strsplit(x,split="T")[[1]][1]), format = "%Y-%m-%d")
dat1$Time.start<-apply(as.data.frame(dat1[,8]),1,FUN=function(x)strsplit(x,split="T")[[1]][2])
dat1$Time.end<-apply(as.data.frame(dat1[,9]),1,FUN=function(x)strsplit(x,split="T")[[1]][2])

#Re-order & remove columns + rename columns:
dat1<-dat1[,c(1,20,10:14,2:7,21,22)]
colnames(dat1)<-c("ID","date","habitat","Total_mass","Subsample_mass","TL.mm","FishNb",
                  "Lat.start","Long.start","Lat.end","Long.end","Depth.start","Depth.end",
                  "Time.start","Time.end")

#Select pelagic herring:
dat1<-dat1[which(dat1[,"habitat"]=="Pelagic"),]

#Calculate fresh weight (FW) of herring numbers based on length-weight relationship 
#(only if total mass of sample was not weighted)
#-->FW (g)=0.315*e^(0.0283*length(mm))
dat1$FW_g<-NA
for(i in 1:nrow(dat1)){
  if(is.na(dat1[i,"Total_mass"])==FALSE){
    dat1[i,"FW_g"]<-NA
  }else{
    dat1[i,"FW_g"]<-(0.315*exp(0.0283*dat1[i,"TL.mm"]))*dat1[i,"FishNb"]
    if(is.na(dat1[i,"TL.mm"])==TRUE)print(paste("Error in line",i,sep=" "))}}
#6 obs are excess=somewhere >1 Total mass was measured or typing error
#a<-dat1[,c("ID","date","Total_mass")]
#a$id<-paste(dat1$ID,dat1$date,sep=" ")
#errors<-as.data.frame(matrix(NA,nrow=0,ncol=4,dimnames=list(NULL,c("ID","date","Total_mass","id"))))
#for(i in 1:length(unique(a$id))){
  #mic<-a[which(a$id==unique(a$id)[i]),]
  #if(length(unique(mic$Total_mass))>1){
   # errors<-rbind(errors,mic)}}
#correct 6 typing errors:
dat1$id<-paste(dat1$ID,dat1$date,sep=" ")
dat1[which(dat1[,"id"]=="Sylt-1 trawl 3 2008-04-07"),"Total_mass"]<-400
dat1[which(dat1[,"id"]=="Sylt-6 trawl 3 2010-05-19"),"Total_mass"]<-112
dat1[which(dat1[,"id"]=="Sylt-9 trawl 1 2014-05-08"),"Total_mass"]<-7897
dat1[which(dat1[,"id"]=="Sylt-6 trawl 4 2015-06-16"),"Total_mass"]<-1223
dat1[which(dat1[,"id"]=="Sylt-1 trawl 5 2016-06-07"),"Total_mass"]<-2541
dat1[which(dat1[,"id"]=="Sylt-2 trawl 4 2019-05-08"),"Total_mass"]<-1011
#Exclude 1 obs: subsample mass indicated but not the total mass (Sylt-2 trawl 6, 2007-06-04)
dat1<-dat1[-which(dat1$id=="Sylt-2 trawl 6 2007-06-04"),]

#Obtain the total fish weight per sample:
totFW<-summaryBy(formula=FW_g~ID+date,dat1,FUN=sum)
dat1<-unique(dat1[,c("ID","date","habitat","Total_mass","Lat.start","Long.start","Lat.end",
                    "Long.end","Depth.start","Depth.end","Time.start","Time.end")])
dat1<-merge(totFW,dat1,by=c("ID","date"),all=TRUE)
dat1$Herring<-NA
for(i in 1:nrow(dat1)){
  if(is.na(dat1[i,"FW_g.sum"])==TRUE){
    dat1[i,"Herring"]<-dat1[i,"Total_mass"]
  }else{
    dat1[i,"Herring"]<-dat1[i,"FW_g.sum"]}}
#Add station column:
dat1$station<-apply(as.data.frame(dat1[,"ID"]),1,FUN=function(x)strsplit(x,split=" ")[[1]][1])
#dat1$station<-sub("ylt-","",dat1$station,fixed=TRUE)
#-->7 stations (?)
#Add year column:
dat1$year<-apply(as.data.frame(dat1[,"date"]),1,FUN=function(x)as.numeric(strsplit(x,split="-")[[1]][1]))
#Add month column:
dat1$month<-apply(as.data.frame(dat1[,"date"]),1,FUN=function(x)as.numeric(strsplit(x,split="-")[[1]][2]))

#Convert fish biomass data into biomass per unit of volume###
#EXCEL-formula to derive trawling distance [in m] per trawl based on start+end coordinates:
#=6371*ARCCOS(COS(BOGENMASS(90-C3))*COS(BOGENMASS(90-E3))+SIN(BOGENMASS(90-C3))*SIN(BOGENMASS(90-E3))*COS(BOGENMASS(D3-F3)))*1000
#-->BOGENMASSfunction~transforms angle in ° into angle in radius: 2pi/360°=x/angle in °
#=>x= angle in ° *2pi/360
#Formula converted to R:
#6371*acos(cos((90-dat1$Lat.start)*2*pi/360)*cos((90-dat1$Lat.end)*2*pi/360)+sin((90-dat1$Lat.start)*2*pi/360)*sin((90-dat1$Lat.end)*2*pi/360)*cos((dat1$Long.start-dat1$Long.end)*2*pi/360))*1000
dist<-function(LAT.START,LAT.END,LONG.START,LONG.END){
  distance<-6371*acos(cos((90-LAT.START)*2*pi/360)*cos((90-LAT.END)*2*pi/360)+sin((90-LAT.START)*2*pi/360)*sin((90-LAT.END)*2*pi/360)*cos((LONG.START-LONG.END)*2*pi/360))*1000
  return(distance)}
#Function test based on calculated excel-data
#dist(55.05555,55.0512,8.438433,8.4534)-1068.998#OK
#dist(54.940833,54.948167,8.399833,8.416333)-1332.498#OK
#dist(55.036,55.0316,8.463283,8.46875)-600.6206#OK
#Trawl opening: 7m*3m
#Sampled Volume=distance*7*3
dat1$samplingDist<-dist(dat1$Lat.start,dat1$Lat.end,dat1$Long.start,dat1$Long.end)
#-->some sampling distances are HUGE: >4.6km
#dat1[which(dat1[,"samplingDist"]>4600000),"date"]#4 sampling dates: "2016-12-14" "2016-02-22" "2016-04-05" "2016-06-07"

#Calculate the sampling volume [in m3] per trawl:
dat1$samplingVol<-dat1$samplingDist*7*3

#Derive fish biomass per unit of volume [in mg/m3]:
dat1$herring.mg_m3<-dat1$Herring*1000/dat1$samplingVol

#Convert fish biomass into biomass carbon [mg C/m3] based on method of De la Vega et al. (2018):
#ash-free DW/FW=0.17 & c/ash-free DW=0.58 (ratios were initially published in Remmert, 1978)
dat1$herring.carbon<-dat1$herring.mg_m3*0.17*0.58

#Add sampling dates at which herring was absent:
#Import df with sampling metadata (2007-2020):
alldates<-fread(file = "Fish/SamplingDates_metadata.csv", na.strings = "", dec = "," , data.table = FALSE)
alldates$Date<-as.Date(alldates$Date, format="%d.%m.%Y")
#Select metadata for pelagic hauls & time period 2007-2019:
alldates<-alldates[which(alldates[,"number of hauls (pelagic)"]==1),]
alldates<-alldates[which(alldates[,"Year"]<=2019),]
alldates<-alldates[-which(is.na(alldates[,"Date"])==TRUE),]#remove some random NA value in date column

#Merge metadata file with fish data:
colnames(alldates)[1:4]<-c("year","month","date","station")#rename columns
dat<-merge(alldates[,c("date","station")],dat1[,c("year","month","date","station","herring.carbon")],
           by=c("date","station"),all=TRUE)
#Detect some typing error:
#-->sampling dates with station of dat1 that were not mentioned in metadata
missingdates<-c()
for(i in 1:nrow(dat1)){
  samplingdate<-paste(dat1[i,"date"], dat1[i,"station"],sep=" ")
  if(length(which(paste(alldates[,"date"],alldates[,"station"],sep=" ")==samplingdate))==0){
    print(dat1[i,c("date","station","herring.carbon")])
    missingdates<-c(missingdates,i)}}

#-->add this missing date to metadata df
alldates<-alldates[,1:4]
alldates<-rbind(alldates,dat1[missingdates,c("year","month","date","station")])
#Retry merging metadata & fish data:
dat<-merge(alldates[,c("year","month","date","station")],dat1[,c("date","station","herring.carbon")],
           by=c("date","station"),all=TRUE)
#Replace NA values of herring.carbon with absence (=0)
dat[which(is.na(dat[,"herring.carbon"])==TRUE),"herring.carbon"]<-0

###################################################################################
###CALCULATE PARAMETER ESTIMATES
source("R functions/get_estimates_ONEparam.R")
stations.her<-c("Sylt-2","Sylt-4","Sylt-6","Sylt-9","Sylt-10")
#-->remove st 1 & 8=outside the bight
data.herring<-get_estimates_ONEparam(dat=dat,timeseries=2007:2019,
                       stations=stations.her,colname.param="herring.carbon")[[2]]

#Difference in annual herring carbon stocks WITH/WITHOUT absence observations (includes data from st1 & st8):
#[1] 0.69552590 0.25748741 0.09546993 0.04143811 0.70335850 1.19660014 0.99048638 0.32917995
#[9] 0.05670073 3.35595741 1.58272001 0.12967984 0.25656713

###################################################################################
####PLOT & EXPLORE DATA
df<-data.frame(year=rep(data.herring$year,2),season=c(data.herring$season,rep("annual",nrow(data.herring))),
               Herring=c(data.herring$seasonal_AVG,data.herring$annual_AVG))

ggplot(df,aes(x=year,y=Herring,color=season,shape=season))+
  theme_classic()+
  geom_point(size=3)+
  geom_smooth(mapping=aes(y=Herring,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  scale_shape_manual("", values=c(15,16,17,18,19),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  ylab("mg C/m3")+
  ggtitle("Clupea harengus (2007-2019)")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="slategray3"))

###################################################################################
####################################################################################
#SEASONALITY OF HERRING IN SRB: 16.12.2008-18.01.2011
setwd(plot.working_dir)
HER2009_11<-dat[min(which(dat$date=="2008-12-16")):max(which(dat$date=="2011-01-18")),]
HER2009_11<-get_estimates_ONEparam(dat=HER2009_11,timeseries=2008:2011,
                                     stations=stations.her,colname.param="herring.carbon")[[1]]
HER2009_11<-HER2009_11[-which(apply(HER2009_11[3:7],1,sum)==0),]
HER2009_11$Date<-paste(HER2009_11$year,HER2009_11$month,sep="-")
HER2009_11$Date<-factor(HER2009_11$Date,levels=HER2009_11$Date)

pdf(file="Herring.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(HER2009_11,aes(x=Date,y=monthly_AVG))+
  theme_minimal()+
  #geom_point(color="darkcyan")+
  geom_line(aes(group=1),size=0.5,color="darkcyan")+
  ylab("mg C/m^3")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"))
# close the graphical device:
dev.off() 
##Smelt (2007-2019)#######################################################################
#dat2<-fread(file = "Fish/data_smelt.csv", na.strings = "", dec = "," , data.table = FALSE)

##Sprat (2007-2019)
#dat3<-fread(file = "Fish/data_sprat.csv", na.strings = "", dec = "," , data.table = FALSE)
