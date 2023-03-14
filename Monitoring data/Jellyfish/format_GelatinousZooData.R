###########################################################
#ANALYSE & FORMAT GELATINOUS ZOOPLANKTON TIME SERIES OF SRB
###########################################################
#-->sampled at List Reede (2009-2021)

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
datall.list<-list()#list to store all imported data files
datall.names<-vector()#vector to store data file names
timeseries<-c(2009:2021)
for(i in 1:length(timeseries)){
  dataname<-paste("Gelatinous_zoo_Reede",timeseries[i],sep="")  
  filename<-paste("Langzeitdaten_J.Rick/Gelatinous zooplankton/",dataname,".csv",sep="")
  dat<-fread(file = filename, na.strings = "", dec = "," , data.table = FALSE)
  if(i<9)colnames(dat)[c(2,3)]<-c("Date","Water_Depth")#rename certain columns
  if(i>=9)colnames(dat)[c(2,5)]<-c("Date","Water_Depth")
  dat<-dat[,-1]#remove unnecessary columns
  dat$year<-timeseries[i]
  dat$station<-"Reede"
  colnames(dat)<-sub(".1","",colnames(dat),fixed=TRUE)
  datall.colnames<-c(datall.colnames,colnames(dat))
  datall.names<-c(datall.names,dataname)
  datall.list[[i]]<-dat}
names(datall.list)<-datall.names

#Merge all data files
datall.colnames<-unique(datall.colnames)#unique column names of all data files
colnames<-datall.colnames[c(20,21,1,27,28,2,22,35,3:19,25,26,29,30,23,24,31:34)]#reorder column names
datall<-as.data.frame(matrix(NA,nrow=10000,ncol=length(colnames),dimnames=list(NULL,colnames)))
nrow.df<-0
for(i in 1:length(datall.list)){
  mydata<-datall.list[[i]]
  for(u in 1:length(colnames)){
    if(length(which(colnames(mydata)==colnames[u]))>0){
      datall[c((nrow.df+1):(nrow.df+nrow(mydata))),colnames[u]]<-mydata[,colnames[u]]}}
  nrow.df<-nrow.df+nrow(mydata)}
datall<-datall[c(1:nrow.df),1:29]#remove empty rows and columns

#rename spp columns for R
allspp<-strsplit(colnames(datall)[c(9:29)],split=" ")
spp.names<-vector()
for(i in 1:length(allspp)){
  myspp<-allspp[[i]]
  myspp<-sub(",","_",myspp,fixed=TRUE)
  myspp<-sub("sp.","_sp",myspp,fixed=TRUE)
  myname<-paste(myspp[1],myspp[2],sep="")
  spp.names<-c(spp.names,myname)}
colnames(datall)[9:29]<-spp.names

#Make date column understandable for R
datall$Date<- as.Date (datall$Date , format = "%d.%m.%Y")

#Add month column
datall$month<-NA
for(i in 1:nrow(datall))datall[i,"month"]<-strsplit(as.character(datall[i,"Date"]),split="-")[[1]][2]
datall$month<-as.numeric(datall$month)

#Convert Nb/3 net hauls into Nb/m^3
#Import df with sampled volumes
Vols<-fread(file = "Langzeitdaten_J.Rick/Gelatinous zooplankton/SampledVolume.csv", na.strings = "", dec = "," , data.table = FALSE)
Vols<-Vols[-which(Vols[,"SampledVolume.m3"]==0),]#remove dates where no sampling took place
Vols$date<-apply(as.data.frame(Vols[,"Date"]),1,FUN=function(x)paste(strsplit(x,split=".",fixed=TRUE)[[1]][3],strsplit(x,split=".",fixed=TRUE)[[1]][2],
                                                                     strsplit(x,split=".",fixed=TRUE)[[1]][1],sep="-"))
for(i in 1:nrow(datall)){
  myVol<-Vols[which(Vols[,"date"]==as.character(datall[i,"Date"])),"SampledVolume.m3"]
  if(length(myVol)==0)myVol<-37.92 #if the date is not indicated in volume-df, use the standard volume
  abundances<-apply(as.data.frame(datall[i,9:29]),2,FUN=function(x)x/myVol)
  datall[i,9:29]<-unname(abundances)}
           
####################################################################################
###SPECIES LIST AND SPECIES GROUPING
#-->Produce a list of all sampled species to facilitate the selection
#of those species that will be taken for the Food web model

#Import list of all sampled gelatinous zooplankton species 
spp.GelZOO<-fread(file = "Langzeitdaten_J.Rick/Gelatinous zooplankton/Gelatinous_zoo_SPPlist.csv", na.strings = "", dec = "," , data.table = FALSE)
spp.GelZOO<-unique(spp.GelZOO)
for(i in 1:nrow(spp.GelZOO)){
  spp<-strsplit(spp.GelZOO[i,"species"],split=" ")[[1]]
  col<-strsplit(spp.GelZOO[i,"column"],split=" ")[[1]]
  if(col[2]=="sp."){
    spp.GelZOO[i,"species"]<-spp[1]
  }else{spp.GelZOO[i,"species"]<-paste(spp[1],spp[2],sep=" ")}
  for(u in 1:length(col)){
    if(u==1){spp.GelZOO[i,"column"]<-col[u]
    }else{spp.GelZOO[i,"column"]<-paste(spp.GelZOO[i,"column"],col[u],sep="")}}
  spp.GelZOO[i,"column"]<-sub("sp.","_sp",spp.GelZOO[i,"column"],fixed=TRUE)
  spp.GelZOO[i,"column"]<-sub(",","_",spp.GelZOO[i,"column"],fixed=TRUE)}
#-->21 species

#Reorder species list according to column order of "data.ZOO"
a<-data.frame(column=spp.names,id=c(1:length(spp.names)))
spp.GelZOO<-merge(a,spp.GelZOO,by="column")
spp.GelZOO<-spp.GelZOO[order(spp.GelZOO[,"id"],decreasing=FALSE),]

#Count the number of missing data points per species (=NA values)
spp.GelZOO$NAs<-apply(datall[,c(9:29)],2,FUN=function(x)length(which(is.na(x)==TRUE)))

#Export spp list:
#write_xlsx(spp.GelZOO, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/Langzeitdaten_J.Rick/Gelatinous zooplankton/GelZooplanktonSpecies.xlsx")

####################################################################################
###SPECIES GROUPING

#Import df with grouping variable
groups.GelZOO<-fread(file = "Langzeitdaten_J.Rick/Gelatinous zooplankton/GROUP.csv", na.strings = "", dec = "," , data.table = FALSE)

#Re-order group df according to column order of datall
#Note=this step seems to be superflu...
ids<-data.frame(column=colnames(datall)[9:(8+length(spp.names))],id=c(1:length(spp.names)))
groups.GelZOO<-merge(groups.GelZOO,ids)
groups.GelZOO<-groups.GelZOO[order(groups.GelZOO$id,decreasing=FALSE),]

#Excluded species/taxonomic groups:
#-->rare or seldomly recorded species
OUTs<-which(groups.GelZOO[,"Group"]=="OUT")+8#9 species/taxonomic groups

#comb jellies:
cteno<-which(groups.GelZOO[,"Group"]=="comb jelly")+8 #5 species
#-->exclude low abundant 
#hydrozoans:
hydro<-which(groups.GelZOO[,"Group"]=="hydrozoan")+8 #10 species

#Preliminary selection:
datall2<-cbind(data.frame(year=datall$year, month=datall$month, date=datall$Date),
               datall[,c(cteno,hydro)])
spp.names<-colnames(datall)[c(cteno,hydro)]                    

###################################################################################
###CALCULATE PARAMETER ESTIMATES
source("R functions/get_estimates_MULTIparam.R")
data.GelZOO<-get_estimates_MULTIparam(dat=datall,timeseries=2009:2021,
                                          stations=c("ST1"),colnames.params=spp.names)[[2]]



###################################################################################
###EXPLORE DATA

#Explore abundances of the selected zooplankton groups in SRB: which are the most dominant?
#mean.GelZoos<-apply(data.GelZOO[,c(3:17)],2,mean,na.rm=TRUE)#mean annual average per gelatinous zooplankton group
#df<-data.frame(species=spp.names,mean=unname(mean.GelZoos))
#df<-df[order(df[,"mean"],decreasing=FALSE),]
#ggplot(df,aes(x=species,y=mean))+
#  theme_minimal()+
#  ylab("Abundance\n[Nb/m^3]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=9,color="slategray3"))

#Explore seasonal patterns of gelatinous zooplankton abundances + temporal trends:
#response<-paste(spp.names,"seasonal_AVG",sep="_")
#response<-response[-c(11,13)]#remove M. haeckelii & Leuckartiara sp. from the series of plots

#Comb jellies
#df1<-data.GelZOO[,c("year","season",response[1:5])]
#df1<-data.frame(year=rep(df1$year,length(cteno)),season=rep(df1$season,length(cteno)),
#                         species=c(rep("Mnemiopsis leidyi",nrow(df1)), rep("Beroe cucumis",nrow(df1)),
#                                   rep("Pleurobrachia pileus",nrow(df1)),rep("Bolinopsis sp.",nrow(df1)),
#                                   rep("juvenile ctenophores",nrow(df1))),
#                abundance=c(df1[,3],df1[,6],df1[,4],df1[,7],df1[,5]))
#df1$species<-factor(df1$species,levels=c("Mnemiopsis leidyi","Beroe cucumis","Pleurobrachia pileus",
#                                         "juvenile ctenophores","Bolinopsis sp."))
#ggplot(df1,aes(x=year,y=abundance,color=season,shape=season))+
#    theme_classic()+
#    geom_point(size=3)+
#    geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df1,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#    scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
#                       breaks=c("winter","spring","summer","autumn"))+
#    scale_shape_manual("", values=c(15,16,17,18),
#                       breaks=c("winter","spring","summer","autumn"))+
#    ylab("Abundance\n[Nb/m^3]")+
#    theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#          axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#          axis.text.y=element_text(face="bold",color="slategray3"),
#          axis.title.x=element_blank(),
#          axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species, scales="free")
  
#Hydrozoans (except M. haeckelii & Leuckartiara sp.as they are very rare)
#df2<-data.GelZOO[,c("year","season",response[6:length(response)])]
#df2<-data.frame(year=rep(df2$year,(length(hydro)-2)),season=rep(df2$season,(length(hydro)-2)),
#                species=c(rep("Sarsia tubulosa",nrow(df2)), rep("Phialella sp.",nrow(df2)),
#                          rep("Rathkea octopuntata",nrow(df2)),
#                          rep("Eucheilota sp.",nrow(df2)),rep("Clytia sp.",nrow(df2)),
#                          rep("Bougainvillia sp.",nrow(df2)),rep("Nemopsis sp.",nrow(df2)),
#                          rep("Obelia geniculata",nrow(df2))),
#                abundance=c(df2[,3],df2[,7],df2[,8],df2[,10],df2[,6],df2[,5],df2[,9],df2[,4]))
#df2$species<-factor(df2$species,levels=c("Rathkea octopuntata","Sarsia tubulosa","Phialella sp.","Eucheilota sp.",
#                                         "Clytia sp.","Bougainvillia sp.","Nemopsis sp.",
#                                         "Obelia geniculata"))
#ggplot(df2,aes(x=year,y=abundance,color=season,shape=season))+
#  theme_classic()+
#  geom_point(size=3)+
#  geom_smooth(mapping=aes(y=abundance,x=year,color=season),data=df2,method="loess",size=0.5, linetype="dashed",se=FALSE)+
#  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod"),
#                     breaks=c("winter","spring","summer","autumn"))+
#  scale_shape_manual("", values=c(15,16,17,18),
#                     breaks=c("winter","spring","summer","autumn"))+
#  ylab("Abundance\n[Nb/m^3]")+
#  theme(#plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#    axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
#    axis.text.y=element_text(face="bold",color="slategray3"),
#    axis.title.x=element_blank(),
#    axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  facet_wrap(~species, scales="free")

############################################################################################
###FINAL SELECTION
hydro<-hydro[-c(which(hydro==20),which(hydro==23))] #remove M. haeckelii & Leuckartiara sp.from hydromedusea group as they are very rare
cteno<-cteno[-c(which(cteno==13),which(cteno==17))]#remove juvenile ctenophores & Bolinopsis sp.

##############################################################################################
###CONVERT ABUNDANCE [Nb/m^3] INTO CARBON BIOMASS [mg C/m^3] BASED ON LENGTH-WEIGHT RELATIONS
#-->for all hydrozoan species for which no average size is available, it is assumed an average size of the compartment group
assumedL.hydro<-mean(c(9,4.5,5.5,6.405))

#Import df with LW-relationships:
LWR<-fread(file = "Langzeitdaten_J.Rick/Gelatinous zooplankton/LWR.csv", na.strings = "", dec = "," , data.table = FALSE)                    

#convert abundances of each species into biomass carbon
CBiom<-datall[,c("year","station","month","Date")]
colnames(CBiom)[4]<-"date"
CBiom<-cbind(CBiom,datall[,c(hydro,cteno)])
for(i in 1:nrow(LWR)){
  myspp<-LWR[i,"Column"]
  LWR.myspp<-parse(text=LWR[i,"LWR"])
  L<-LWR[i,"meanL"]
  if(is.na(L)==TRUE)L<-assumedL.hydro
  C.myspp<-eval(LWR.myspp)
  print(C.myspp)
  C.ratio<-LWR[i,"C/W-ratio"] #(if necessary) adjust biomass estimate based on DW/WW and C/WW ratios
  if(is.na(C.ratio)==TRUE)C.ratio<-1
  print(C.ratio)
  data.myspp<-which(colnames(CBiom)==myspp)
  CBiom[,data.myspp]<-CBiom[,data.myspp]*C.myspp*C.ratio}


###################################################################################
###CALCULATE PARAMETER ESTIMATES
datall2<-data.frame(year=datall$year, month=datall$month, date=datall$Date,
                    M.leidyi=CBiom[,"M.leidyi"],
                    P.pileus=CBiom[,"P.pileus"],
                    B.cucumis=CBiom[,"B.cucumis"],
                    Hydromedusae=apply(CBiom[,colnames(datall)[hydro]],1,sum,na.rm=TRUE))
spp.names<-c("M.leidyi","P.pileus","B.cucumis","Hydromedusae")

data.GelZOO<-get_estimates_MULTIparam(dat=datall2,timeseries=2009:2021,
                                      stations=c("ST1"),colnames.params=spp.names)[[2]]

#Explore abundances of the selected zooplankton groups in SRB: which are the most dominant?
#mean.GelZoos<-apply(data.GelZOO[,c(3:6)],2,mean,na.rm=TRUE)#mean annual average per gelatinous zooplankton group
#df<-data.frame(species=spp.names,mean=unname(mean.GelZoos))
#df<-df[order(df[,"mean"],decreasing=FALSE),]
#ggplot(df,aes(x=species,y=mean))+
#  theme_minimal()+
#  ylab("Abundance\n[mg C/m^3]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=9,color="slategray3"))
#
#Explore seasonal patterns of gelatinous zooplankton abundances + temporal trends:
response<-paste(spp.names,"seasonal_AVG",sep="_")
title<-c("Mnemiopsis leidyi","Pleurobrachia pileus","Beroe cucumis","Hydromedusae")
for(i in 1:length(spp.names)){
  df<-data.GelZOO[,c("year","season",response[i])]
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

###################################################################################
####################################################################################
#SEASONALITY OF GEL. ZOOPLANKTON IN SRB: 11.05.2009-21.02.2011
setwd(plot.working_dir)
GelZOO2009_11<-CBiom
GelZOO2009_11$Hydromedusae<-apply(GelZOO2009_11[,5:12],1,sum)
GelZOO2009_11<-GelZOO2009_11[which(GelZOO2009_11[,"date"]<="2011-02-21"),]
gelzoo.names<-colnames(GelZOO2009_11)[13:16]
GelZOO2009_11<-get_estimates_MULTIparam(dat=GelZOO2009_11,timeseries=2009:2011,
                                     stations=c("ST1"),
                                     colnames.params=gelzoo.names)[[1]]
GelZOO2009_11<-GelZOO2009_11[-which(GelZOO2009_11[,"Nobs_ST1"]==0),]
GelZOO2009_11$Date<-paste(GelZOO2009_11$year,GelZOO2009_11$month,sep="-")
GelZOO2009_11$Date<-factor(GelZOO2009_11$Date,levels=GelZOO2009_11$Date)

#Plot seasonality
gelzoo2009_11<-as.data.frame(matrix(NA,nrow=0,ncol=3,
                                 dimnames=list(NULL,c("Time","Comp","Carbon"))))
 
for(i in 1:length(gelzoo.names)){
  name<-paste(gelzoo.names[i],"monthly_AVG",sep="_")
  df<-data.frame(Time=GelZOO2009_11$Date,
                 Comp=rep(gelzoo.names[i],nrow(GelZOO2009_11)),
                 Carbon=GelZOO2009_11[,which(colnames(GelZOO2009_11)==name)])
  gelzoo2009_11<-rbind(gelzoo2009_11,df)}   
gelzoo2009_11$Time<-factor(gelzoo2009_11$Time,levels=GelZOO2009_11$Date)
gelzoo2009_11$Comp<-factor(gelzoo2009_11$Comp,levels=c("M.leidyi","P.pileus","B.cucumis","Hydromedusae"))
#Store data:
write.csv(gelzoo2009_11,file=paste(plot.working_dir,"Biomass/jellyfish.csv",sep=""))

cols<-c("red","goldenrod","chocolate","peru","orange",
        "slateblue3","mediumseagreen","plum","violetred","cyan")
pdf(file="Jellyfish.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(gelzoo2009_11,aes(x=Time,y=Carbon,group=Comp,color=Comp))+
  theme_minimal()+
  geom_point()+
  geom_line(size=0.5)+
  ylab("mg C/m^3")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"),
        legend.title = element_blank())
# close the graphical device:
dev.off() 
#########################################################################
####CARBON CONTRIBUTION OF EACH COMPARTMENT TO TOTAL CARBON CARBON
##############ACROSS SEASONS 2009/2010 (winter 2008/09-winter 2010/11)
#-->NA values mean that compartment seasonal biomass = 0 (dividing something by 0 is not possible)
setwd(working_dir)
source("R functions/carbon_contributions.R")

spp.names<-vector(length=length(5:12))
a<-1
for(i in 5:12){
  print(i)
if(length(strsplit(colnames(CBiom)[i],split="_",fixed=TRUE)[[1]])>1){
spp.names[a]<-paste(strsplit(colnames(CBiom)[i],split="_",fixed=TRUE)[[1]][1], "sp.", sep=" ")
}else if(length(strsplit(colnames(CBiom)[i],split="_",fixed=TRUE)[[1]])==1){
  spp.names[a]<-colnames(CBiom)[i]}
  a<-a+1}
CARBcontr.GelZOO<-spptocompC(dat=CBiom,stations="1",
                          species=list(spp.names,c("M.leidyi","P.pileus","B.cucumis")),
                          columns= list(colnames(CBiom)[5:12],c(colnames(CBiom)[13],colnames(CBiom)[14],colnames(CBiom)[15])),
                          comps=c("Hydromedusae","ctenophores"))
#Save and plot results
setwd(plot.working_dir)
write_xlsx(CARBcontr.GelZOO,"Jellyfish.xlsx")
pdf(file="JellyfishSpeciesComposition.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
plot(plot.Ccontributions(CARBcontr.GelZOO))
# close the graphical device:
dev.off()
