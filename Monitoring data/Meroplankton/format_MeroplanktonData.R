###################################################
#ANALYSE & FORMAT MEROPLANKTON TIME SERIES OF SRB
###################################################
#-->sampled at List Ferry Terminal (2007-2021) & Station 1 (2014-2021)

#set working directory
working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/"
plot.working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Community composition/"
setwd(working_dir)
#load packages
library(data.table)
library(doBy)
library(ggplot2)
library(writexl)

##############################################################################
###IMPORT & FORMAT DATA

datall.colnames<-vector()#column names of all data files

#Ferry terminal:
datall_FT<-list()#list to store all imported data files
datall_FT.names<-vector()#vector to store data file names
timeseries<-c(2007:2021)
for(i in 1:length(timeseries)){
dataname<-paste("Meroplankton_FT_",timeseries[i],sep="")  
filename<-paste("Langzeitdaten_J.Rick/Meroplankton/",dataname,".csv",sep="")
dat<-fread(file = filename, na.strings = "N/A", dec = "," , data.table = FALSE)
colnames(dat)[c(2:5)]<-c("Date","Latitude","Longitude","Water_Depth")#rename certain columns
dat<-dat[-1,-c(1,6)]#remove unnecessary rows+columns
dat$year<-timeseries[i]
dat$station<-"FT"
datall.colnames<-c(datall.colnames,colnames(dat))
datall_FT.names<-c(datall_FT.names,dataname)
datall_FT[[i]]<-dat}
names(datall_FT)<-datall_FT.names

#Station 1:
datall_ST1<-list()#list to store all imported data files
datall_ST1.names<-vector()#vector to store data file names
timeseries<-c(2014:2021)
for(i in 1:length(timeseries)){
  dataname<-paste("Meroplankton_ST1_",timeseries[i],sep="")  
  filename<-paste("Langzeitdaten_J.Rick/Meroplankton/",dataname,".csv",sep="")
  dat<-fread(file = filename, na.strings = "N/A", dec = "," , data.table = FALSE)
  colnames(dat)[c(2:5)]<-c("Date","Latitude","Longitude","Water_Depth")#rename certain columns
  dat<-dat[-1,-c(1,6)]#remove unnecessary rows+columns
  dat$year<-timeseries[i]
  dat$station<-"ST1"
  datall.colnames<-c(datall.colnames,colnames(dat))
  datall_ST1.names<-c(datall_ST1.names,dataname)
  datall_ST1[[i]]<-dat}
names(datall_ST1)<-datall_ST1.names

#Merge all data files
datall.colnames<-unique(datall.colnames)#unique column names of all data files
colnames<-datall.colnames[c(30,31,1:29,34:37,38,32,33)]#reorder column names
datall<-as.data.frame(matrix(NA,nrow=10000,ncol=38,dimnames=list(NULL,colnames)))
nrow.df<-0
for(i in 1:length(datall_FT)){#add all data files from Ferry Terminal
mydata<-datall_FT[[i]]
for(u in 1:length(colnames)){
  if(length(which(colnames(mydata)==colnames[u]))>0){
  datall[c((nrow.df+1):(nrow.df+nrow(mydata))),colnames[u]]<-mydata[,colnames[u]]}}
nrow.df<-nrow.df+nrow(mydata)}
for(i in 1:length(datall_ST1)){#add all data files from station 1
  mydata<-datall_ST1[[i]]
  for(u in 1:length(colnames)){
    if(length(which(colnames(mydata)==colnames[u]))>0){
      datall[c((nrow.df+1):(nrow.df+nrow(mydata))),colnames[u]]<-mydata[,colnames[u]]}}
  nrow.df<-nrow.df+nrow(mydata)}
datall<-datall[c(1:nrow.df),c(1:35)]#remove empty rows and columns
datall<-datall[-which(datall$Date==""),]
#Rename spp columns for R
allspp<-strsplit(colnames(datall)[c(7:ncol(datall))],split=" ")
spp.names<-vector()
for(i in 1:length(allspp)){
  myspp<-allspp[[i]]
  if(length(myspp)>1){myname<-paste(myspp[1],myspp[2],sep="_")
  }else{myname<-myspp}
  spp.names<-c(spp.names,myname)}
colnames(datall)[7:ncol(datall)]<-spp.names

#Make date column understandable for R
datall$Date<- as.Date (datall$Date , format = "%d.%m.%Y")

#Add month column
datall$month<-NA
for(i in 1:nrow(datall))datall[i,"month"]<-strsplit(as.character(datall[i,"Date"]),split="-")[[1]][2]

#Format abundance/month columns into numeric
for(i in 7:ncol(datall))datall[,i]<-as.numeric(datall[,i])

#Check sampling frequency at station 1 & FT since 2014
#summaryBy(Nobs.sum~station, 
          #data=data.MERO[which(data.MERO[,"year"]>=2014),],FUN=sum)
#ST1: 3508 sampling dates
#FT: 3508 sampling dates

####################################################################################
###SPECIES LIST
#-->Produce a list of all sampled species to facilitate the selection
#of those species that will be taken for the Food web model

#Import list of all sampled meroplankton species 
spp.MERO<-fread(file = "Langzeitdaten_J.Rick/Meroplankton/Meroplankton_SPPlist.csv", na.strings = "", dec = "," , data.table = FALSE)
spp.MERO<-unique(spp.MERO)
#-->nb of species columns (29) does not fit with number of species in list (31):
######spp list: "Pandalus" & "Pandalus Megalopa" are both named "Pandalus" in the column=THE SAME(?)
######spp list: "Crangon Megalopa" & "Crangon" are both named "Crangon" in the column=THE SAME(?)
#=>PROBABLY 29 species=remove excess rows from species list:
spp.MERO<-spp.MERO[-c(which(spp.MERO[,"species"]=="Pandalus"),which(spp.MERO[,"species"]=="Crangon")),]
for(i in 1:nrow(spp.MERO)){
  spp<-strsplit(spp.MERO[i,"column"],split=" ",fixed=TRUE)[[1]]
  for(u in 1:length(spp)){
    if(u==1){spp.MERO[i,"column"]<-spp[u]
    }else{spp.MERO[i,"column"]<-paste(spp.MERO[i,"column"],spp[u],sep="_")}}}

#Reorder species list according to column order of data frame
a<-data.frame(column=spp.names,id=c(1:length(spp.names)))
spp.MERO<-merge(a,spp.MERO,by="column")
spp.MERO<-spp.MERO[order(spp.MERO[,"id"],decreasing=FALSE),c("id","species","column")]

#Count the number of missing data points per species (=NA values)
spp.MERO$NAs<-apply(datall[,c(7:35)],2,FUN=function(x)length(which(is.na(x)==TRUE)))

#Export spp list:
#write_xlsx(spp.MERO, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/Langzeitdaten_J.Rick/Meroplankton/MeroplanktonSpecies.xlsx")

####################################################################################
###SPECIES GROUPING

#Import df with grouping variable
groups.MERO<-fread(file = "Langzeitdaten_J.Rick/Meroplankton/GROUP.csv", na.strings = "", dec = "," , data.table = FALSE)

#Re-order group df according to column order of datall
ids<-data.frame(column=colnames(datall)[7:(6+length(spp.names))],id=c(1:length(spp.names)))
groups.MERO<-merge(groups.MERO,ids)
groups.MERO<-groups.MERO[order(groups.MERO$id,decreasing=FALSE),]

#Excluded species/taxonomic groups:
#-->rare or seldomly recorded species/taxonomic groups
OUTs<-which(groups.MERO[,"Group"]=="OUT")+6#14 species/taxonomic groups
#Bivalve larvae:
biv<-which(groups.MERO[,"Group"]=="bivalve larvae")+6
#Gastropod larvae:
gastr<-which(groups.MERO[,"Group"]=="gastropod larvae")+6
#Cypris larvae:
cyp<-which(groups.MERO[,"Group"]=="cyprid")+6
#Polychaete larvae:
poly<-which(groups.MERO[,"Group"]=="polychaete larvae")+6#5 species
#Echinoderm larvae:
echi<-which(groups.MERO[,"Group"]=="echinoderm larvae")+6
#Tunicate larvae (=appendicularians):
#-->larvae or adults?
tunicates<-which(groups.MERO[,"Group"]=="tunicate")
#Cladoceran larvae:
#-->probably not used for the model
clado<-which(groups.MERO[,"Group"]=="cladoceran larvae")+6

#################################################################################################
##EXPLORE DATA####################
source("R functions/get_estimates_MULTIparam.R")
###
#Polychaetes: seasonal cycle of polychaete larvae at species-level 
#df1<-datall[,c("year","station","month","Date","Spionidae","Scoloplos","Nereis","Lanice","Magelona")]
#spp.names<-c("Spionidae","Scoloplos","Nereis","Lanice","Magelona")
#df1<-get_estimates_MULTIparam(dat=df1,timeseries=2007:2021,
#                              stations=c("ST1","FT"),colnames.params=spp.names)[[2]]
#df1<-data.frame(year=rep(df1$year,length(spp.names)),season=rep(df1$season,length(spp.names)),
#                species=c(rep("Spionidae",nrow(df1)),rep("Scoloplos",nrow(df1)), 
#                          rep("Nereis",nrow(df1)),rep("Lanice",nrow(df1)),
#                          rep("Magelona",nrow(df1))),
#                abundance=c(df1$Spionidae_seasonal_AVG,df1$Scoloplos_seasonal_AVG, 
#                            df1$Nereis_seasonal_AVG,df1$Lanice_seasonal_AVG,df1$Magelona_seasonal_AVG))
#df1$species<-factor(df1$species,levels=c("Spionidae","Lanice","Nereis","Magelona","Scoloplos"))
#ggplot(df1,aes(x=year,y=abundance,color=season,shape=season))+
#  theme_classic()+
#  geom_point(size=3)+
#  geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df1,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  scale_shape_manual("", values=c(15,16,17,18,3),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  ylab("Abundance [Nb/10L]")+
#  theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#    axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#    axis.text.y=element_text(face="bold",color="slategray3"),
#    axis.title.x=element_blank(),
#    axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species,scales="free",ncol=2)
#Abundance fractions:
#sum(apply(datall[,poly],1,sum,na.rm=TRUE),na.rm=TRUE)#total abundance:61710/10L
#sum(apply(datall[,spp.names[2:length(spp.names)]],1,sum,na.rm=TRUE),na.rm=TRUE)#3122/10L
#-->Scoloplos/Nereis/Lanice/Magelona make up 5.059148% of total polychaete larval abundance
#-->Sabine´s suggestion: Lanice/Nereis are important in terms of carbon
#=>ONLY EXCLUDE SCOLOPLOS FROM POLYCHAETE LARVAE GROUP DUE TO LOW ABUNDANCE+REL.SMALL SIZE:300-900um

###
#Echinoderm larvae: seasonal cycle of bipinnaria & echinopluteus
#df2<-datall[,c("year","station","month","Date","Bipinnaria","Echinopluteus")]
#spp.names<-c("Bipinnaria","Echinopluteus")
#df2<-get_estimates_MULTIparam(dat=df2,timeseries=2007:2021,
#                              stations=c("ST1","FT"),colnames.params=spp.names)[[2]]
#df2<-data.frame(year=rep(df2$year,length(spp.names)),season=rep(df2$season,length(spp.names)),
#                species=c(rep("Bipinnaria",nrow(df2)),rep("Echinopluteus",nrow(df2))),
#                abundance=c(df2$Bipinnaria_seasonal_AVG,df2$Echinopluteus_seasonal_AVG))
#ggplot(df2,aes(x=year,y=abundance,color=season,shape=season))+
#  theme_classic()+
#  geom_point(size=3)+
#  geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df2,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  scale_shape_manual("", values=c(15,16,17,18,3),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  ylab("Abundance [Nb/10L]")+
#  theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#    axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#    axis.text.y=element_text(face="bold",color="slategray3"),
#    axis.title.x=element_blank(),
#    axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species,scales="free",ncol=1)

############################################################################################
###FINAL SELECTION
poly<-c(which(colnames(datall)=="Spionidae"),which(colnames(datall)=="Lanice"),
        which(colnames(datall)=="Nereis"),which(colnames(datall)=="Magelona"))#exclude Scoloplos from polychaete larvae group

##############################################################################################
###CONVERT ABUNDANCE [Nb/10L] INTO CARBON BIOMASS [mg C/m^3] BASED ON LENGTH-WEIGHT RELATIONS

#Import df with LW-relationships:
LWR<-fread(file = "Langzeitdaten_J.Rick/Meroplankton/LWR.csv", na.strings = "", dec = "," , data.table = FALSE)                    
LWR[,"LWR"]<-sub(",",".",LWR[,"LWR"],fixed=TRUE)

#convert abundances of each species into biomass carbon
CBiom<-datall[,c("year","station","month","Date")]
colnames(CBiom)[4]<-"date"
CBiom$gastropods<-apply(datall[,gastr],1,sum,na.rm=TRUE)
CBiom<-cbind(CBiom,datall[,c(poly,biv)])
for(i in 1:nrow(LWR)){
  myspp<-LWR[i,"Column"]
  LWR.myspp<-parse(text=LWR[i,"LWR"])
  L<-LWR[i,"meanL"]
  #if(is.na(L)==TRUE)L<-assumedL.hydro
  C.myspp<-eval(LWR.myspp)
  C.ratio<-LWR[i,"C/W-ratio"] #(if necessary) adjust biomass estimate based on DW/WW and C/WW ratios
  if(is.na(C.ratio)==TRUE)C.ratio<-1
  #print(paste(myspp,C.myspp*C.ratio))
  data.myspp<-which(colnames(CBiom)==myspp)
  if(length(which(LWR[,"Column"]==myspp))>1){
    CBiom[,(ncol(CBiom)+1)]<-CBiom[,data.myspp]*C.myspp*C.ratio
  }else{
  CBiom[,data.myspp]<-CBiom[,data.myspp]*C.myspp*C.ratio}}

CBiom[,"Spionidae"]<-apply(CBiom[,c(11,12)],1,mean,na.rm=TRUE)#calculate Spionidae biomass based on the average of the two applied LWRs

#Convert biomass estimates from ugC/10L into mgC/m3:
#10L=0,01m3 ; 1000ug=1mg
#*100-->/m3 ; *0.001-->mg => conversion factor of 0.1
CBiom[,5:10]<-apply(CBiom[,5:10],2,FUN=function(x)x*0.1)

################################################################################################
###CALCULATE PARAMETER ESTIMATES
datall2<-data.frame(year=datall$year,station=datall$station,month=datall$month, date=datall$Date,
                    bivalves=CBiom[,"Bivalvia"],
                    polychaetes=apply(CBiom[,c("Spionidae","Lanice","Magelona","Nereis")],1,sum,na.rm=TRUE),
                    gastropods=CBiom[,"gastropods"])
spp.names<-c("polychaetes","bivalves","gastropods")

data.MERO<-get_estimates_MULTIparam(dat=datall2,timeseries=2007:2021,
                                stations=c("1","FT"),colnames.params=spp.names)[[2]]



###################################################################################
###EXPLORE DATA

#Explore abundances of the selected zooplankton groups in SRB: which are the most dominant?
#mean.meros<-apply(data.MERO[,c(4:6)],2,mean,na.rm=TRUE)#mean annual average per meroplankton group
#df<-data.frame(species=spp.names,mean=unname(mean.meros))
#df<-df[order(df[,"mean"],decreasing=FALSE),]
#ggplot(df,aes(x=species,y=mean))+
#  theme_minimal()+
#  ylab("Abundance\n[mg C/m^3]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=9,color="slategray3"))


#Explore seasonal patterns of meroplankton abundances + temporal trends:
response<-paste(spp.names,"seasonal_AVG",sep="_")
title<-c("Polychaete larvae","Bivalve larvae","Gastropod larvae")
for(i in 1:length(spp.names)){
  df<-data.MERO[,c("year","season",response[i])]
  colnames(df)<-c("year","season","Abundance")
  plot.title<-title[i]
  myplot<-ggplot(df,aes(x=year,y=Abundance,color=season,shape=season))+
    theme_classic()+
    geom_point(size=3)+
    geom_smooth(mapping=aes(y=Abundance,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
    scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
                       breaks=c("winter","spring","summer","autumn"))+
    scale_shape_manual("", values=c(15,16,17,18),
                       breaks=c("winter","spring","summer","autumn"))+
    ylab("Abundance\n[mg C/m^3]")+
    ggtitle(plot.title)+
    theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
          axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
          axis.text.y=element_text(face="bold",color="slategray3"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(face="bold", size=9,color="slategray3"))
  print(myplot)}

#########################################################################
#########################################################################
#SEASONALITY OF MEROPLANKTON IN SRB: 11.12.2008-23.02.2011
setwd(plot.working_dir)
MERO2009_11<-CBiom[,1:10]
MERO2009_11$Polychaetes<-apply(MERO2009_11[,6:9],1,sum)
MERO2009_11<-MERO2009_11[which(MERO2009_11[,"date"]>="2008-12-11"),]
MERO2009_11<-MERO2009_11[which(MERO2009_11[,"date"]<="2011-02-23"),]
mero.names<-colnames(MERO2009_11)[5:11]
MERO2009_11<-get_estimates_MULTIparam(dat=MERO2009_11,timeseries=2008:2011,
                                     stations=c("FT","ST1"),
                                     colnames.params=mero.names)[[1]]
MERO2009_11<-MERO2009_11[-which(MERO2009_11[,"Nobs_FT"]==0),]
MERO2009_11$Date<-paste(MERO2009_11$year,MERO2009_11$month,sep="-")
MERO2009_11$Date<-factor(MERO2009_11$Date,levels=MERO2009_11$Date)

#Plot seasonality
spp<-colnames(MERO2009_11)[5:10]
groups<-colnames(MERO2009_11)[c(5,10,11)]
mero.sp2009_11<-as.data.frame(matrix(NA,nrow=0,ncol=3,
                                    dimnames=list(NULL,c("Time","Spp","Carbon"))))
mero2009_11<-as.data.frame(matrix(NA,nrow=0,ncol=3,
                                 dimnames=list(NULL,c("Time","Comp","Carbon"))))
for(i in 1:length(spp)){
  name<-strsplit(spp[i],split="_",fixed=TRUE)[[1]][1]
  df<-data.frame(Time=MERO2009_11$Date,
                 Spp=rep(name,nrow(MERO2009_11)),
                 Carbon=MERO2009_11[,which(colnames(MERO2009_11)==spp[i])])
  mero.sp2009_11<-rbind(mero.sp2009_11,df)}  
for(i in 1:length(groups)){
  name<-strsplit(groups[i],split="_",fixed=TRUE)[[1]][1]
  df<-data.frame(Time=MERO2009_11$Date,
                 Comp=rep(name,nrow(MERO2009_11)),
                 Carbon=MERO2009_11[,which(colnames(MERO2009_11)==groups[i])])
  mero2009_11<-rbind(mero2009_11,df)}   

mero.sp2009_11$Time<-factor(mero.sp2009_11$Time,levels=MERO2009_11$Date)
mero2009_11$Time<-factor(mero2009_11$Time,levels=MERO2009_11$Date)
mero2009_11$Comp<-factor(mero2009_11$Comp,levels=c("Bivalvia","Polychaetes","gastropods"))
#Store data:
write.csv(mero2009_11,file=paste(plot.working_dir,"Biomass/meroplankton.csv",sep=""))

#Species-level:
spec<-c("Bivalvia","Spionidae","Lanice","Nereis","Magelona","gastropods")
cols<-c("mediumseagreen","hotpink","plum","violetred","darkmagenta","cyan")
mero.sp2009_11$Spp<-factor(mero.sp2009_11$Spp,levels=spec)
pdf(file="MeroplanktonSpecies.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(mero.sp2009_11,aes(x=Time,y=Carbon,group=Spp,color=Spp))+
  theme_minimal()+
  geom_point()+
  geom_line(size=0.5)+
  ylab("mg C/m^3")+
  scale_color_manual("",values=cols,breaks=spec)+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"),
        legend.title = element_blank())
# close the graphical device:
dev.off() 
#Compartment-level:
pdf(file="Meroplankton.pdf",         # File name
    width = 6, height = 3, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(mero2009_11,aes(x=Time,y=Carbon,group=Comp,color=Comp))+
  theme_bw()+
  #geom_point()+
  geom_line(size=0.75)+
  scale_color_manual("Meroplankton",values=c("lightcoral","lightblue","hotpink4"),
                     breaks=c("Bivalvia","Polychaetes","gastropods"),
                     labels=c("Bivalvia","Polychaeta","Gastropoda"))+
  ylab("mg C/m^3")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_blank(),#element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"))
# close the graphical device:
dev.off() 
#########################################################################
####CARBON CONTRIBUTION OF EACH COMPARTMENT TO TOTAL CARBON CARBON
##############ACROSS SEASONS 2009/2010 (winter 2008/09-winter 2010/11)
#-->NA values mean that compartment seasonal biomass = 0 (dividing something by 0 is not possible)
setwd(working_dir)
source("R functions/carbon_contributions.R")

CARBcontr.MERO<-spptocompC(dat=CBiom,stations=c("FT","ST1"),
                             species=list(colnames(CBiom)[6:9],"Bivalvia","Gastropoda"),
                             columns= list(colnames(CBiom)[6:9],colnames(CBiom)[10],colnames(CBiom)[5]),
                             comps=c("polychaetes","bivalves","gastropods"))
#Save and plot results:
setwd(plot.working_dir)
write_xlsx(CARBcontr.MERO,"Meroplankton.xlsx")
pdf(file="MeroplanktonSpeciesComposition.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
plot(plot.Ccontributions(CARBcontr.MERO))
# close the graphical device:
dev.off() 
