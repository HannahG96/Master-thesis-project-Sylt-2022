#################################################
#ANALYSE & FORMAT ZOOPLANKTON TIME SERIES OF SRB 
#################################################
#-->sampled at

#set working directory
working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/"
plot.working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Community composition/"
setwd(working_dir)

#load packages
library(ggplot2)
library(writexl)
library(data.table)

##############################################################################
###IMPORT & FORMAT DATA

source("R functions/import_pangaea_file.R")#function to convert tab file into data frame
folder<-"Zooplankton_P.Martens/"
a<-rep("ListReede_zooplankton_",5)
timeseries<-c(2007:2011)
filenames<-vector()#Create vector of the filenames to be imported
for(i in 1:length(timeseries)){
  myfilename<-paste(folder,a[i],timeseries[i],".tab",sep="")
  filenames<-c(filenames,myfilename)}

datall.list<-list()#list to store all imported data files
datall.names<-vector()#vector to attribute a name to each data file
datall.colnames<-vector()#column names of all data files
stations<-rep("List Reede",5)
for(i in 1:length(filenames)){
  mydata.name<-strsplit(filenames[i],split="/",fixed=TRUE)[[1]][2]
  mydata.name<-strsplit(mydata.name,split=".",fixed=TRUE)[[1]][1]#Name of the data file
  datall.names<-c(datall.names,mydata.name)
  dat<-import_pangaea_file(myfile=filenames[i])
  dat$year<-timeseries[i]#column that indicates sampling year of the data file
  if(length(which(colnames(dat)=="Event"))==0)dat$Event<-stations[i]#if necessary, column that indicates sampling station of the data file
  datall.colnames<-c(datall.colnames,colnames(dat))
  datall.list[[i]]<-dat}
names(datall.list)<-datall.names

#Merge all data files
datall.colnames<-unique(datall.colnames)#unique column names of all data files
colnames<-datall.colnames[c(55,54,1:53)]#reorder column names
datall<-as.data.frame(matrix(NA,nrow=10000,ncol=length(colnames),dimnames=list(NULL,colnames)))
nrow.df<-0
for(i in 1:length(datall.list)){
  mydata<-datall.list[[i]]
  for(u in 1:length(colnames)){
    if(length(which(colnames(mydata)==colnames[u]))>0){
      datall[c((nrow.df+1):(nrow.df+nrow(mydata))),colnames[u]]<-mydata[,colnames[u]]}}
  nrow.df<-nrow.df+nrow(mydata)}
datall<-datall[c(1:nrow.df),]#remove empty rows

#Make date column understandable for R
datall$Date_Time<- as.Date (datall$Date_Time , format = "%Y-%m-%d")

#Add month column
#Add month+season column
datall$month<-NA
for(i in 1:nrow(datall))datall[i,"month"]<-strsplit(as.character(datall[i,"Date_Time"]),split="-")[[1]][2]

#Format abundance columns into numeric
for(i in 5:ncol(datall))datall[,i]<-as.numeric(datall[,i])
####Note: in zooplankton tab file 2011 -->missing values, NAs introduced

####################################################################################
###SPECIES LIST 
#-->Produce a list of all sampled species to facilitate the selection
#of those species that will be taken for the Food web model

#Extract all species names of the sampled zooplankton from tab files
source("R functions/get_species.R")
spp.ZOO<-as.data.frame(matrix(NA,nrow=0,ncol=2,dimnames=list(NULL,c("species","column"))))
for(i in 1:length(filenames)){
  spp<-get_species(myfile=filenames[i])
  spp.ZOO<-rbind(spp.ZOO,spp)}
spp.ZOO<-unique(spp.ZOO)
spp.ZOO$genus<-NA
for(i in 1:nrow(spp.ZOO))spp.ZOO[i,"genus"]<-strsplit(spp.ZOO[i,"species"],split=" ")[[1]][1]
spp.ZOO$genus<-sub(",","",spp.ZOO$genus,fixed=TRUE)
#-->51 species ; 34 genera

#Reorder species list according to column order of data frame
spp.names<-colnames(datall)[5:55]
a<-data.frame(column=spp.names,id=c(1:length(spp.names)))
spp.ZOO<-merge(a,spp.ZOO,by="column")
spp.ZOO<-spp.ZOO[order(spp.ZOO[,"id"],decreasing=FALSE),c("id","species","genus","column")]

#Count the number of missing data points per species (=NA values)
spp.ZOO$NAs<-apply(datall[,c(5:55)],2,FUN=function(x)length(which(is.na(x)==TRUE)))

#Export spp list:
#write_xlsx(spp.ZOO, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/Zooplankton_P.Martens/ZooplanktonSpecies.xlsx")

####################################################################################
###SPECIES GROUPING

#Import df with grouping variable
groups.ZOO<-fread(file = "Zooplankton_P.Martens/GROUP.csv", na.strings = "", dec = "," , data.table = FALSE)

#Re-order group df according to column order of datall
#Note=this step seems to be superflu...
ids<-data.frame(column=colnames(datall)[5:(4+length(spp.names))],id=c(1:length(spp.names)))
groups.ZOO<-merge(groups.ZOO,ids)
groups.ZOO<-groups.ZOO[order(groups.ZOO$id,decreasing=FALSE),]

#Identify data columns corresponding to each functional group

#Excluded species/taxonomic groups:
#-->overlaps with meroplankton, microplankton & gel. zooplankton data sets, rare or seldomly recorded species
OUTs<-which(groups.ZOO[,"Group"]=="OUT")+4#22 species/taxonomic groups

#Copepod nauplii
naups<-which(groups.ZOO[,"Group"]=="copepod nauplii")+4#5 species/taxonomic groups

#Calanoid copepods:
#OPTION1:copepodites+adults (-->option rejected after data exploration=go with the simplest method)
#cals1<-c(which(groups.ZOO[,"Group"]=="calanoid copepod"),
         #which(groups.ZOO[,"Group"]=="calanoid copepodites"))+4#15 species/taxonomic groups
#OPTION2:only adults
cals<-which(groups.ZOO[,"Group"]=="calanoid copepod")+4#10 species/taxonomic groups

#Harpacticoid copepods
harp<-which(groups.ZOO[,"Group"]=="harpacticoid copepod")+4#1 taxonomic group

#Oithona similis (cyclopoid copepod)
#OPTION1:adults+copepodites (-->option rejected, treat cyclopoids like calanoids)
#Oith1<-c(which(groups.ZOO[,"Group"]=="cyclopoid copepod"),
        #which(groups.ZOO[,"Group"]=="cyclopoid copepodites"))+4
#OPTION2:only adults
Oith<-which(groups.ZOO[,"Group"]=="cyclopoid copepod")+4

#cladocerans
clado<-which(groups.ZOO[,"Group"]=="cladoceran")+4

#rotifers
rot<-which(groups.ZOO[,"Group"]=="rotifer")+4

#tunicates
tun<-which(groups.ZOO[,"Group"]=="tunicate")+4


###EXPLORE DATA##############################################################################
source("R functions/get_estimates_MULTIparam.R")
####
#calanoid copepods at species-level
#df1<-data.frame(year=datall$year,month=datall$month,date=datall$Date_Time,
                #calanoids=apply(datall[,cals],1,sum,na.rm=TRUE),
#                Acartia_sp=apply(datall[,c("Acartia_sp._f","Acartia_sp._m")],1,sum,na.rm=TRUE),
#                T.longicornis=apply(datall[,c("T._longicornis_f","T._longicornis_m")],1,sum,na.rm=TRUE),
#                C.hamatus=apply(datall[,c("C._hamatus_f","C._hamatus_m")],1,sum,na.rm=TRUE),
#                P.parvus=apply(datall[,c("P._parvus_f","P._parvus_m")],1,sum,na.rm=TRUE),
#                P.elongatus=apply(datall[,c("P._elongatus_f","P._elongatus_m")],1,sum,na.rm=TRUE))
spp.names<-c(#"calanoids",
  "Acartia_sp","T.longicornis","C.hamatus","P.parvus","P.elongatus")
#df1<-get_estimates_MULTIparam(dat=df1,timeseries=2007:2011,
#                              stations="ST1",colnames.params=spp.names)[[2]]
#df1<-data.frame(year=rep(df1$year,length(spp.names)),season=rep(df1$season,length(spp.names)),
#                species=c(#rep("Calanoid copepods",nrow(df1)),
#                  rep("Acartia sp.",nrow(df1)),rep("Temora longicornis",nrow(df1)), 
#                          rep("Centropages hamatus",nrow(df1)),rep("Paracalanus parvus",nrow(df1)),
#                          rep("Pseudocalanus elongatus",nrow(df1))),
#                abundance=c(#df1$calanoids_seasonal_AVG,
#                  df1$Acartia_sp_seasonal_AVG,df1$T.longicornis_seasonal_AVG, 
#                            df1$C.hamatus_seasonal_AVG,df1$P.parvus_seasonal_AVG,df1$P.elongatus_seasonal_AVG))
#df1$species<-factor(df1$species,levels=c(#"Calanoid copepods",
#  "Acartia sp.","Temora longicornis","Centropages hamatus","Paracalanus parvus","Pseudocalanus elongatus"))
#ggplot(df1,aes(x=year,y=abundance,color=season,shape=season))+
#  theme_classic()+
#  geom_point(size=3)+
#  geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df1,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  scale_shape_manual("", values=c(15,16,17,18,3),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  ylab("Abundance [Nb/m^3]")+
#  theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#    axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#    axis.text.y=element_text(face="bold",color="slategray3"),
#    axis.title.x=element_blank(),
#    axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species,scales="free",ncol=2)

######
#cladocerans at species-level
#df2<-data.frame(year=datall$year,month=datall$month,date=datall$Date_Time,
               #cladocerans=apply(datall[,clado],1,sum,na.rm=TRUE),
#               Podon_sp=datall[,c("Podon_sp.")],
#               Evadne_sp=datall[,c("Evadne_sp.")])
#spp.names<-c(#"cladocerans",
#             "Podon_sp","Evadne_sp")
#df2<-get_estimates_MULTIparam(dat=df2,timeseries=2007:2011,
#                             stations="ST1",colnames.params=spp.names)[[2]]
#df2<-data.frame(year=rep(df2$year,length(spp.names)),season=rep(df2$season,length(spp.names)),
#               species=c(#rep("Cladocerans",nrow(df2)),
#                 rep("Podon sp.",nrow(df2)),rep("Evadne sp.",nrow(df2))),
#               abundance=c(#df2$cladocerans_seasonal_AVG,
#                           df2$Podon_sp_seasonal_AVG,
#                           df2$Evadne_sp_seasonal_AVG))
#df2$species<-factor(df2$species,levels=c(#"Cladocerans",
#  "Podon sp.","Evadne sp."))
#ggplot(df2,aes(x=year,y=abundance,color=season,shape=season))+
#  theme_classic()+
#  geom_point(size=3)+
#  geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df2,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  scale_shape_manual("", values=c(15,16,17,18,3),
#                     breaks=c("winter","spring","summer","autumn","annual"))+
#  ylab("Abundance [Nb/m^3]")+
#  theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#    axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#    axis.text.y=element_text(face="bold",color="slategray3"),
#   axis.title.x=element_blank(),
#    axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species,scales="free",ncol = 1)

##############################################################################################
###CONVERT ABUNDANCE [Nb/m^3] INTO CARBON BIOMASS [mg C/m^3] BASED ON LENGTH-WEIGHT RELATIONS
#-->for all copepod species for which no average C/W-ratio is available, it is assumed an average C/W-ratio of the compartment group
assumedDW.copepods<-mean(c(0.465,0.497,0.433,0.453))

CBiom<-data.frame(year=datall$year, month=datall$month, date=datall$Date_Time,
                  tunicates=datall[,tun],
                  Oithona_sp=apply(datall[,Oith],1,sum,na.rm=TRUE),
                  harpacticoids=datall[,harp],
                  Acartia_sp=apply(datall[,c("Acartia_sp._f","Acartia_sp._m")],1,sum,na.rm=TRUE),
                  T.longicornis=apply(datall[,c("T._longicornis_f","T._longicornis_m")],1,sum,na.rm=TRUE),
                  C.hamatus=apply(datall[,c("C._hamatus_f","C._hamatus_m")],1,sum,na.rm=TRUE),
                  P.parvus=apply(datall[,c("P._parvus_f","P._parvus_m")],1,sum,na.rm=TRUE),
                  P.elongatus=apply(datall[,c("P._elongatus_f","P._elongatus_m")],1,sum,na.rm=TRUE),
                  Podon_sp=datall[,c("Podon_sp.")],
                  Evadne_sp=datall[,c("Evadne_sp.")])
copepods<-colnames(CBiom)[5:11]
#Import df with LW-relationships
LWR<-fread(file = "Zooplankton_P.Martens/LWR.csv", na.strings = "", dec = "," , data.table = FALSE)                    
LWR$LWR<-sub(",",".",LWR$LWR,fixed=TRUE)

#convert abundances of each species into biomass carbon
for(i in 1:nrow(LWR)){
  myspp<-LWR[i,"Column"]
  LWR.myspp<-parse(text=LWR[i,"LWR"])
  L<-LWR[i,"meanL"]
  C.myspp<-eval(LWR.myspp)
  #print(C.myspp)
  C.ratio<-LWR[i,"C/W-ratio"] #(if necessary) adjust biomass estimate based on C/W-ratios
  if(length(which(copepods==myspp))!=0 &&
     is.na(C.ratio)==TRUE){ #DW/WW of copepods is assumed 15%
    C.ratio<-0.15*assumedDW.copepods
  }else if(length(which(copepods==myspp))){
    C.ratio<-C.ratio*0.15
  }else if(is.na(C.ratio)==TRUE){#convert carbon biomass of tunicates+cladocerans from ugC to mgC
    C.ratio<-0.001}
  data.myspp<-which(colnames(CBiom)==myspp)
  CBiom[,data.myspp]<-CBiom[,data.myspp]*C.myspp*C.ratio}
write_xlsx(CBiom,path="C:/Hannah/Zooplankton2007-2012.xlsx")


###############################################################################################
###CALCULATE PARAMETER ESTIMATES
datall2<-data.frame(year=datall$year, month=datall$month, date=datall$Date_Time,
                    copepods=apply(CBiom[,c("Acartia_sp","C.hamatus",
                                             "T.longicornis","P.parvus",
                                             "P.elongatus","Oithona_sp","harpacticoids")],1,sum,na.rm=TRUE),
                    cladocerans=apply(CBiom[,c("Podon_sp","Evadne_sp")],1,sum,na.rm=TRUE),
                    tunicates=CBiom[,"tunicates"])
spp.names<-c("copepods", "cladocerans", "tunicates")

data.ZOO<-get_estimates_MULTIparam(dat=datall2,timeseries=2007:2011,
                                stations="ST1",colnames.params=spp.names)[[2]]
#write_xlsx(data.ZOO,path="C:/Hannah/Zooplankton2007-2012_seasonalAVG.xlsx")

###################################################################################
###EXPLORE DATA

#Explore abundances of the selected zooplankton groups in SRB: which are the most dominant?
#mean.zoos<-apply(data.ZOO[,c(3:5)],2,mean,na.rm=TRUE)#mean annual average per zooplankton group
#df<-data.frame(species=spp.names,mean=unname(mean.zoos))
#df<-df[order(df[,"mean"],decreasing=FALSE),]
#ggplot(df,aes(x=species,y=mean))+
#  theme_minimal()+
#  ylab("Abundance\n[mg C/m^3]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=9,color="slategray3"))


#Explore seasonal patterns of zooplankton abundances + temporal trends:
response<-paste(spp.names,"seasonal_AVG",sep="_")
title<-c("Copepods",
         #"Copepod nauplii",
         "Cladocerans",
         #"Rotifers",
         "Tunicates")
for(i in 1:length(spp.names)){
  df<-data.ZOO[,c("year","season",response[i])]
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


####################################################################################
#SEASONALITY OF ZOOPLANKTON IN SRB: 01.12.2008-21.02.2011

ZOO2009_11<-CBiom
ZOO2009_11$Copepods<-apply(ZOO2009_11[,5:11],1,sum)
ZOO2009_11$Cladocerans<-apply(ZOO2009_11[,12:13],1,sum)
ZOO2009_11<-ZOO2009_11[which(ZOO2009_11[,"date"]>="2008-12-01"),]
ZOO2009_11<-ZOO2009_11[which(ZOO2009_11[,"date"]<="2011-02-21"),]
zoo.names<-colnames(ZOO2009_11)[4:15]
ZOO2009_11<-get_estimates_MULTIparam(dat=ZOO2009_11,timeseries=2008:2011,
                                     stations=c("ST1"),
                                     colnames.params=zoo.names)[[1]]
ZOO2009_11<-ZOO2009_11[-which(ZOO2009_11[,"Nobs_ST1"]==0),]
ZOO2009_11$Date<-paste(ZOO2009_11$year,ZOO2009_11$month,sep="-")
ZOO2009_11$Date<-factor(ZOO2009_11$Date,levels=ZOO2009_11$Date)

#Plot seasonality
spp<-colnames(ZOO2009_11)[4:13]
groups<-colnames(ZOO2009_11)[c(4,14,15)]
zoo.sp2009_11<-as.data.frame(matrix(NA,nrow=0,ncol=3,
                  dimnames=list(NULL,c("Time","Spp","Carbon"))))
zoo2009_11<-as.data.frame(matrix(NA,nrow=0,ncol=3,
                                    dimnames=list(NULL,c("Time","Spp","Carbon"))))
for(i in 1:length(spp)){
  name<-strsplit(spp[i],split="_",fixed=TRUE)[[1]][1]
  df<-data.frame(Time=ZOO2009_11$Date,
                 Spp=rep(name,nrow(ZOO2009_11)),
                 Carbon=ZOO2009_11[,which(colnames(ZOO2009_11)==spp[i])])
  zoo.sp2009_11<-rbind(zoo.sp2009_11,df)}  
for(i in 1:length(groups)){
  name<-strsplit(groups[i],split="_",fixed=TRUE)[[1]][1]
  df<-data.frame(Time=ZOO2009_11$Date,
                 Comp=rep(name,nrow(ZOO2009_11)),
                 Carbon=ZOO2009_11[,which(colnames(ZOO2009_11)==groups[i])])
  zoo2009_11<-rbind(zoo2009_11,df)}   
  
zoo.sp2009_11$Time<-factor(zoo.sp2009_11$Time,levels=ZOO2009_11$Date)
zoo2009_11$Time<-factor(zoo2009_11$Time,levels=ZOO2009_11$Date)
zoo2009_11$Comp<-factor(zoo2009_11$Comp,levels=c("Copepods","Cladocerans","tunicates"))
#Store data:
write.csv(zoo2009_11,file=paste(plot.working_dir,"Biomass/zooplankton.csv",sep=""))

#Species-level:
spec<-c("Acartia","T.longicornis","C.hamatus","P.elongatus","P.parvus",
        "Oithona","harpacticoids","Podon","Evadne","tunicates")
cols<-c("red","goldenrod","chocolate","peru","orange",
        "slateblue3","mediumseagreen","plum","violetred","cyan")
zoo.sp2009_11$Spp<-factor(zoo.sp2009_11$Spp,levels=spec)

#Plot Results:
setwd(plot.working_dir)
pdf(file="MesozooplanktonSpecies.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(zoo.sp2009_11,aes(x=Time,y=Carbon,group=Spp,color=Spp))+
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
pdf(file="Mesozooplankton.pdf",         # File name
    width = 6, height = 3, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(zoo2009_11,aes(x=Time,y=Carbon,group=Comp,color=Comp))+
  theme_bw()+
  #geom_point()+
  geom_line(size=0.75)+
  scale_color_manual("Mesozooplankton",values=c("gold","firebrick1","lightslateblue"),
                     breaks=c("Copepods","Cladocerans","tunicates"),
                     labels=c("copepods","cladocerans","tunicates"))+
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

CARBcontr.ZOO<-spptocompC(dat=CBiom,stations="1",
    species=list(colnames(CBiom)[5:11],colnames(CBiom)[12:13],"O.dioica"),
    columns= list(colnames(CBiom)[5:11],colnames(CBiom)[12:13],colnames(CBiom)[4]),
    comps=c("copepods","cladocerans","tunicates"))
#Save and plot results
setwd(plot.working_dir)
write_xlsx(CARBcontr.ZOO,"Mesozooplankton.xlsx")
pdf(file="MesozooplanktonSpeciesComposition.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
plot(plot.Ccontributions(CARBcontr.ZOO))
# close the graphical device:
dev.off() 
