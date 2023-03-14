###################################################
#ANALYSE & FORMAT HYDROCHEMISTRY TIME SERIES OF SRB
###################################################

#set working directory
#set working directory
working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/"
plot.working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Community composition/"
setwd(working_dir)
#load packages
library(data.table)
library(doBy)
library(ggplot2)
library(writexl)

#######################################################################################
###IMPORT & FORMAT HYDROCHEMISTRY DATA

#dat<-fread(file = "Langzeitdaten_J.Rick/Hydrochemistry/SRB_Hydrochemistry_Timeseries.csv", na.strings = "", dec = "," , data.table = FALSE)
#head(dat)
dat1<-fread(file = "Langzeitdaten_J.Rick/Hydrochemistry/SRB_Hydrochemistry_Timeseries.csv", na.strings = "", dec = "," , data.table = FALSE)
dat2<-fread(file = "Langzeitdaten_J.Rick/Hydrochemistry/SRB_Hydrochemistry_2020.csv", na.strings = "", dec = "," , data.table = FALSE)
dat3<-fread(file = "Langzeitdaten_J.Rick/Hydrochemistry/SRB_Hydrochemistry_2021.csv", na.strings = "", dec = "," , data.table = FALSE)

colnames(dat2)<-colnames(dat3)<-c("Date","OUT","station","salinity","SST","pH","PO4",
                                  "Si","NH4","NO2","NO3","SPM","chlorophyll a")
dat2$latitude<-dat2$longitude<-NA
dat3$latitude<-dat3$longitude<-NA
dat2[which(dat2[,"station"]=="FÃ¤hre"),"station"]<-"Ferry Terminal"
dat3[which(dat3[,"station"]=="FÃ¤hre"),"station"]<-"Ferry Terminal"
dat2<-dat2[,-which(colnames(dat2)=="OUT")]
dat3<-dat3[,-which(colnames(dat3)=="OUT")]
dat2<-dat2[,colnames(dat1)]
dat3<-dat3[,colnames(dat1)]
dat<-rbind(dat1,dat2,dat3)

#Make date column understandable for R
dat$Date<- as.Date (dat$Date , format = "%d.%m.%Y")
#Measurement stations:
unique(dat$station)#5 stations
#Add year & month column
dat$year<-dat$month<-NA
for(i in 1:nrow(dat)){
dat[i,"year"]<-strsplit(as.character(dat[i,"Date"]),split="-")[[1]][1]
dat[i,"month"]<-strsplit(as.character(dat[i,"Date"]),split="-")[[1]][2]}
dat$year<-as.numeric(dat$year)
dat$month<-as.numeric(dat$month)

#Check sampling frequency at station 1, 4 & FT since 2009
#summaryBy(Nobs.sum~station, 
          #data=data.SST[which(data.SST[,"year"]>=2009),],FUN=sum)
#ST1: 3552 sampling dates
#ST4: 1169 sampling dates
#FT: 482 sampling dates

#Select data set for timeperiod 1980-2021 of station 1 & 4
#dat<-dat[which(dat[,"year"]>=1980),]
#dat<-dat[-c(which(dat[,"station"]=="Ferry Terminal"),
           # which(dat[,"station"]=="2"),which(dat[,"station"]=="3")),]
#data start: 08.01.1990; data end: 29.11.2021

##################################################################################
##SST
###################

#check for NA values
nrow(dat[which(is.na(dat[,"SST"])==TRUE),c("Date","station","year")])#34 obs

#Calculate (i)annual AVG + MAX + MIN & (ii)seasonal AVG + SD + %SD of SRB
source("R functions/get_estimates_ONEparam.R")
data.SST<-get_estimates_ONEparam(dat=dat,timeseries=1985:2021,
                                 stations=c("Ferry Terminal","1","4"),colname.param="SST")[[2]]

#Plot the data:
ggplot(data.SST,aes(x=year,y=seasonal_AVG,color=season))+
  theme_classic()+
  geom_line()+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
                     breaks=c("winter","spring","summer","autumn"))+
  geom_line(mapping=aes(y=annual_AVG,x=year),color="black",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MAX,x=year),color="plum",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MIN,x=year),color="plum",linetype="dashed")+
  ylab("SST (°C)")+
  theme(axis.text.x=element_text(face="bold",angle=45,color="slategray3",size=8),
        axis.text.y=element_text(face="bold",color="slategray3",size=7),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_text(face="bold"))
#no temperature sampling in 1983-->remove winter average:
data.SST[13,"seasonal_AVG"]<-NA
#+no temperature sampling in winter 1999

#Investigate global warming trends
#1.annual average
ModelSST.1<-lm(annual_AVG~year,data=data.SST)
summary(ModelSST.1)#-->significant with Rsquared of 0.2725
#2.spring average
ModelSST.2<-lm(seasonal_AVG~year,data=data.SST[which(data.SST[,"season"]=="spring"),])
summary(ModelSST.2)#-->significant with Rsquared of 0.1206
#3.summer average
ModelSST.3<-lm(seasonal_AVG~year,data=data.SST[which(data.SST[,"season"]=="summer"),])
summary(ModelSST.3)#-->significant with Rsquared of 0.2751
#4.autumn average
ModelSST.4<-lm(seasonal_AVG~year,data=data.SST[which(data.SST[,"season"]=="autumn"),])
summary(ModelSST.4)#-->significant with Rsquared of 0.08402
#5.winter average
ModelSST.5<-lm(seasonal_AVG~year,data=data.SST[which(data.SST[,"season"]=="winter"),])
summary(ModelSST.5)#-->significant with Rsquared of 0.1217

#Plot temperature anomaly based average annual T° during 1985-2021
annualA<-unique(data.SST$annual_AVG)-mean(data.SST$annual_AVG,na.rm=TRUE)
springA<-data.SST[which(data.SST[,"season"]=="spring"),"seasonal_AVG"]-
  mean(data.SST[which(data.SST[,"season"]=="spring"),"seasonal_AVG"],na.rm=TRUE)
summerA<-data.SST[which(data.SST[,"season"]=="summer"),"seasonal_AVG"]-
  mean(data.SST[which(data.SST[,"season"]=="summer"),"seasonal_AVG"],na.rm=TRUE)
autumnA<-data.SST[which(data.SST[,"season"]=="autumn"),"seasonal_AVG"]-
  mean(data.SST[which(data.SST[,"season"]=="autumn"),"seasonal_AVG"],na.rm=TRUE)
winterA<-c(data.SST[which(data.SST[,"season"]=="winter"),"seasonal_AVG"]-
  mean(data.SST[which(data.SST[,"season"]=="winter"),"seasonal_AVG"],na.rm=TRUE))
anomalies<-data.frame(year=rep(1985:2021,5),
                      season=c(rep("annual",length(1985:2021)),rep("spring",length(1985:2021)),
                               rep("summer",length(1985:2021)),rep("autumn",length(1985:2021)),
                               rep("winter",length(1985:2021))),
                      anomaly=c(annualA,springA,summerA,autumnA,winterA))
anomalies$year<-factor(anomalies$year,levels=as.character(1985:2021))
negs<-pos<-anomalies
negs[which(negs[,"anomaly"]>0),"anomaly"]<-0
pos[which(pos[,"anomaly"]<0),"anomaly"]<-0
negs$season<-factor(negs$season, levels=c("annual","winter","spring","summer","autumn"))
pos$season<-factor(pos$season, levels=c("annual","winter","spring","summer","autumn"))
anomalies$season<-factor(anomalies$season, levels=c("annual","winter","spring","summer","autumn"))
ggplot(anomalies, aes(x = year, y = anomaly)) + 
  theme_bw()+
  geom_bar(data = negs,position="dodge", stat = "identity",fill="darkcyan") +
  geom_bar(data = pos,position="dodge", stat = "identity",fill="gold") + 
  #geom_errorbarh(color="gray48",position =position_dodge(.9), height=0.25)+
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.3)+
  ylab("Temperature anomaly [°C]")+
  xlab("Year")+
  theme(panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line( size=.1, color="grey" ),
        legend.title = element_blank(),
        axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x = element_text(face="bold", size=10, color="grey28"),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~season,ncol=1)

#Create df of seasonal SSTs for the period winter 2008/9 - winter 2010/11:
seasons<-c("2009 winter","2009 spring","2009 summer","2009 autumn","2010 winter",
           "2010 spring","2010 summer","2010 autumn","2011 winter")
SST_2009_11<-data.SST[c(which(data.SST$year==2009),which(data.SST$year==2010),
                        which(data.SST$year==2011)),c("year","season","seasonal_AVG",
                                                      "seasonal_SD","seasonal_pctSD")]
SST_2009_11$season<-paste(SST_2009_11$year,SST_2009_11$season,sep=" ")
for(i in 1:length(seasons)){ #chronologically order data frame:
  if(i==1){
    SST2009_11<-SST_2009_11[which(SST_2009_11$season==seasons[i]),c("season","seasonal_AVG",
                                                                 "seasonal_SD","seasonal_pctSD")]
  }else{
    SST2009_11<-rbind(SST2009_11,SST_2009_11[which(SST_2009_11$season==seasons[i]),c("season","seasonal_AVG",
                                                                "seasonal_SD","seasonal_pctSD")])}}
#Export data of seasonal SSTs as Excel file:
#write_xlsx(SST2009_11, 
     # path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/SST.xlsx")

#####################################################################################
##SALINITY
###################
setwd(plot.working_dir)
#Get monthly average Salinity values (01.12.2008-25.02.2011):
data.SAL<-dat[which(dat$Date<="2011-02-25"),]
data.SAL<-data.SAL[which(data.SAL$Date>="2008-12-01"),]
#Plot temporal trends of salinity:
data.SAL$time<-as.character(data.SAL$Date)
data.SAL$time<-factor(data.SAL$time,levels=unique(data.SAL$time))

ggplot(data.SAL,aes(x=time,y=salinity))+
  theme_minimal()+
  geom_point(color="gray48")+
  geom_line(aes(group=1),size=0.5,color="gray48")+
  ylab("PSU")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=90,size=4,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"))

data.SAL<-get_estimates_ONEparam(dat=data.SAL,timeseries=2008:2011,
                                 stations=c("Ferry Terminal","1","4"),colname.param="salinity")[[1]]
data.SAL<-data.SAL[-which(is.na(data.SAL$monthly_AVG)==TRUE),]
data.SAL$season<-paste(data.SAL$year,data.SAL$month,sep="-")
data.SAL$season<-factor(data.SAL$season,levels=data.SAL$season)
#Plot the data:
pdf(file="Salinity.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(data.SAL,aes(x=season,y=monthly_AVG))+
  theme_minimal()+
  geom_point(color="gray48")+
  geom_line(aes(group=1),size=0.5,color="gray48")+
  ylab("PSU")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"))
# close the graphical device:
dev.off() 

data.SAL<-get_estimates_ONEparam(dat=dat,timeseries=1985:2021,
                                 stations=c("Ferry Terminal","1","4"),colname.param="salinity")[[2]]

#Plot the data:
ggplot(data.SAL,aes(x=year,y=seasonal_AVG,color=season))+
  theme_classic()+
  geom_line()+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
                     breaks=c("winter","spring","summer","autumn"))+
  geom_line(mapping=aes(y=annual_AVG,x=year),color="black",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MAX,x=year),color="plum",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MIN,x=year),color="plum",linetype="dashed")+
  ylab("Salinity")+
  theme(axis.text.x=element_text(face="bold",angle=45,color="slategray3",size=8),
        axis.text.y=element_text(face="bold",color="slategray3",size=7),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_text(face="bold"))
#####################################################################################
##CHL A
###################

colnames(dat)[13]<-"Chla"#rename chlorophyll a column

#check for NA values
dat[which(is.na(dat[,"Chla"])==TRUE),c("Date","station","year")]#4 obs

#Check sampling frequency at station 1, 4 & FT since 2009
#summaryBy(Nobs.sum~station, 
          #data=data.Chla[which(data.Chla[,"year"]>=2009),],FUN=sum)
#ST1: 3552 sampling dates
#ST4: 1165 sampling dates
#FT: 503 sampling dates


#Calculate (i)annual AVG + MAX + MIN & (ii)seasonal AVG + SD + %SD of SRB
data.Chla<-get_estimates_ONEparam(dat=dat,timeseries=2000:2019,
                                  stations=c("1","4"),colname.param="Chla")[[2]]



#Plot the data
ggplot(data.Chla,aes(x=year,y=seasonal_AVG,color=season))+
  theme_classic()+
  geom_line()+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
                     breaks=c("winter","spring","summer","autumn"))+
  geom_line(mapping=aes(y=annual_AVG,x=year),color="black",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MAX,x=year),color="plum",linetype="dashed")+
  #geom_line(mapping=aes(y=annual_MIN,x=year),color="plum",linetype="dashed")+
  ylab("Chlorophyll a [ug/L]")+
  theme(axis.text.x=element_text(face="bold",angle=45,color="slategray3",size=8),
        axis.text.y=element_text(face="bold",color="slategray3",size=7),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_text(face="bold"))

df<-data.frame(year=rep(data.Chla$year,2),
               season=c(data.Chla$season,rep("annual",nrow(data.Chla))),
               param=c(data.Chla$seasonal_AVG,data.Chla$annual_AVG))
ggplot(df,aes(x=year,y=param,color=season,shape=season))+
  theme_classic()+
  geom_point(size=3)+
  geom_smooth(mapping=aes(y=param,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  scale_shape_manual("", values=c(15,16,17,18,3),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  ylab("Chl a [ug/L]")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="slategray3"))
#-->spring outlier in 2010 (?)

#####################################################################################
##SPM=suspended particulate matter
##################################

#check for NA values
dat[which(is.na(dat[,"SPM"])==TRUE),c("Date","station","year")]#9 obs

#Check sampling frequency at station 1, 4 & FT since 2009
#summaryBy(Nobs.sum~station, 
          #data=data.SPM[which(data.SST[,"year"]>=2009),],FUN=sum)
#ST1: 3536 sampling dates
#ST4: 1161 sampling dates
#FT: 503 sampling dates

#Calculate (i)annual AVG + MAX + MIN & (ii)seasonal AVG + SD + %SD of SRB
data.SPM<-get_estimates_ONEparam(dat=dat,timeseries=2007:2019,
                                  stations=c("1","4"),colname.param="SPM")[[2]]

#Plot the data:
ggplot(data.SPM,aes(x=year,y=seasonal_AVG,color=season))+
  theme_classic()+
  geom_line()+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
                     breaks=c("winter","spring","summer","autumn"))+
  geom_line(mapping=aes(y=annual_AVG,x=year),color="black",linetype="dashed")+
 # geom_line(mapping=aes(y=annual_MAX,x=year),color="plum",linetype="dashed")+
 # geom_line(mapping=aes(y=annual_MIN,x=year),color="plum",linetype="dashed")+
  ylab("SPM [mg/L]")+
  theme(axis.text.x=element_text(face="bold",angle=45,color="slategray3",size=8),
        axis.text.y=element_text(face="bold",color="slategray3",size=7),
        axis.title.x=element_text(face="bold"),
        axis.title.y=element_text(face="bold"))
df<-data.frame(year=rep(data.SPM$year,2),
               season=c(data.SPM$season,rep("annual",nrow(data.SPM))),
               param=c(data.SPM$seasonal_AVG,data.SPM$annual_AVG))
ggplot(df,aes(x=year,y=param,color=season,shape=season))+
  theme_classic()+
  geom_point(size=3)+
  geom_smooth(mapping=aes(y=param,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  scale_shape_manual("", values=c(15,16,17,18,3),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  ylab("SPM [mg/L]")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="slategray3"))
#-->winter outlier in 2013 (?)