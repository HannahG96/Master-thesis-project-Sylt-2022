### DETECT THE 10 MOST IMPORTANT CARBON FLOWS IN EACH SEASON ###

#Set working directory:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Results/")

#Load packages:
library(data.table)
library(writexl)
#Import data on seasonal food web flows:
seasonfolders<-c("Spring 2009","Summer 2009","Autumn 2009","Spring 2010","Summer 2010","Autumn 2010")
seasons<-c("2009 spring","2009 summer","2009 autumn",
           "2010 spring","2010 summer","2010 autumn")
for(i in 1:length(seasonfolders)){
  if(i==1){
    dat<-fread(file = paste(seasonfolders[i],"FLOWS_CI95.csv",sep="/"), na.strings = "", dec = "," , 
               data.table = FALSE)
    dat[,"season"]<-seasons[i]
  }else{
      subdat<-fread(file = paste(seasonfolders[i],"FLOWS_CI95.csv",sep="/"), na.strings = "", dec = "," , 
                    data.table = FALSE)
      subdat$season<-seasons[i]
      dat<-rbind(dat,subdat)}}
#Select the 10 largest carbon flows of each season:
BIGS<-list()
for(i in 1:length(seasons)){
  mydat<-dat[which(dat[,"season"]==seasons[i]),]
  MIN<-mydat[order(mydat[,"Q2.5"],decreasing=TRUE),c("Flow","Mean","Q2.5","Q97.5")]
  MIN<-MIN[1:10,]
  MAX<-mydat[order(mydat[,"Q97.5"],decreasing=TRUE),c("Flow","Mean","Q2.5","Q97.5")]
  MAX<-MAX[1:10,]
  bigs<-merge(MIN,MAX,by=c("Flow","Mean","Q2.5","Q97.5"),all=TRUE)
  bigs<-bigs[order(bigs[,"Mean"],decreasing=TRUE),]
  BIGS[[i]]<-bigs}
names(BIGS)<-seasons
write_xlsx(BIGS,"bigFlows.xlsx")
