###########################################################
#ANALYSE & FORMAT PHYTO-/MICROPLANKTON TIME SERIES OF SRB (species-level)
###########################################################

#NOTE: the summed carbon values changed a little bit since I added the command to order
#"group.PhytoMicro"-rows according to "datall"-columns

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

source("R functions/import_pangaea_file.R")#function to convert tab file into data frame
folder<-"Phyto-_Microplankton.spp_J.Rick/"
a<-c(rep("List_phytoplank_biomassC_",4),rep("List_Reede_phytoplank_microzooplank_biomassC_",3),
     rep("List_Ferry_Terminal_phytoplank_microzooplank_biomassC_",3))
timeseries<-c(2007:2013,2011:2013)
filenames<-vector()#Create vector of the filenames to be imported
for(i in 1:length(timeseries)){
  myfilename<-paste(folder,a[i],timeseries[i],".tab",sep="")
  filenames<-c(filenames,myfilename)}

datall.list<-list()#list to store all imported data files
datall.names<-vector()#vector to attribute a name to each data file
datall.colnames<-vector()#column names of all data files
stations<-c(rep(NA,4),rep("List_Reede",3),rep("List_Ferry_Terminal",3))#sampling stations
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
spp.names_columns<-datall.colnames
datall.colnames<-unique(datall.colnames)#unique column names of all data files
colnames<-datall.colnames[c(1,151,2:149,152:228,150)]#reorder column names
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
datall$month<-NA
for(i in 1:nrow(datall))datall[i,"month"]<-strsplit(as.character(datall[i,"Date_Time"]),split="-")[[1]][2]

#Format biomass carbon+month columns into numeric
for(i in 4:ncol(datall))datall[,i]<-as.numeric(datall[,i])

#Select the data from station 1:
datall<-datall[which(datall[,"Event"]=="List_Reede"),]

####################################################################################
###SPECIES LIST 
#-->Produce a list of all sampled species to facilitate the selection
#of those species that will be taken for the Food web model

#Extract all species names of the sampled phyto-/microplankton from tab files
source("R functions/get_species.R")
spp.PhytoMicro<-as.data.frame(matrix(NA,nrow=0,ncol=2,dimnames=list(NULL,c("species","column"))))
for(i in 1:length(filenames)){
spp<-get_species(myfile=filenames[i])
spp.PhytoMicro<-rbind(spp.PhytoMicro,spp)}
spp.PhytoMicro<-unique(spp.PhytoMicro)
spp.PhytoMicro$genus<-NA
for(i in 1:nrow(spp.PhytoMicro))spp.PhytoMicro[i,"genus"]<-strsplit(spp.PhytoMicro[i,"species"],split=" ")[[1]][1]
#-->224 species ; 133 genera

#Reorder species list according to column order of data frame
spp.names<-colnames(datall)[5:228]
a<-data.frame(column=spp.names,id=c(1:length(spp.names)))
spp.PhytoMicro<-merge(a,spp.PhytoMicro,by="column")
spp.PhytoMicro<-spp.PhytoMicro[order(spp.PhytoMicro[,"id"],decreasing=FALSE),c("id","species","genus","column")]

#Count the number of missing data points per species (=NA values)
spp.PhytoMicro$NAs<-apply(datall[,c(5:228)],2,FUN=function(x)length(which(is.na(x)==TRUE)))

#Export spp list:
#write_xlsx(spp.PhytoMicro, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/Phyto-_Microplankton.spp_J.Rick/MicroplanktonSpecies.xlsx")


####################################################################################
###SPECIES GROUPING
#-->Form species groups & pool corresponding data

#PRE-SELECTION##
#Import df with grouping variable
groups.PhytoMicro<-fread(file = "Phyto-_Microplankton.spp_J.Rick/GROUP.csv", na.strings = "", dec = "," , data.table = FALSE)

#Re-order group df according to column order of datall
ids<-data.frame(column=colnames(datall)[5:(4+length(spp.names))],id=c(1:length(spp.names)))
groups.PhytoMicro<-merge(groups.PhytoMicro,ids)
groups.PhytoMicro<-groups.PhytoMicro[order(groups.PhytoMicro$id,decreasing=FALSE),]

#Identify data columns corresponding to each functional group

#Excluded species/taxonomic groups:
#-->rare or seldomly recorded
#eg. amoeba, chlorophytes, some haptophytes, a raphidophyte...
OUTs<-which(groups.PhytoMicro[,"group"]=="OUT")+4#22 species/taxonomic groups

#diatoms:
diatoms<-which(groups.PhytoMicro[,"group"]=="diatom")+4#115 species/taxonomic groups

#dinoflagellates:
#-->dominated by N. scintillans(=green form with photosynthetic symbiont Pedinomonas noctiluca is autotrophic, red form is heterotrophic)
#-->ASK MARTHE WHICH FORM DOMINATES IN SRB!!
dinos<-which(groups.PhytoMicro[,"group"]=="dinoflagellate")+4#56 species/taxonomic groups

#Phaeocystis sp. (haptophyte):
Phaeo<-which(groups.PhytoMicro[,"group"]=="Phaeocystis sp.")+4

#Cryptophyceae:
Crypto<-which(groups.PhytoMicro[,"group"]=="Cryptophyceae")+4
sum(datall[,Crypto],na.rm=TRUE)#1442.35 ug C/L

#Dictyocha speculum:
Dspec<-which(groups.PhytoMicro[,"group"]=="Dictyocha speculum")+4
sum(datall[,Dspec],na.rm=TRUE)#65.425 ug C/L

#heterotrophic flagellates:
#-->choanoflagellates & E. tripartita are mainly heterotroph
#-->unclear whether undetermined flagellates can be classified as heterotrophic
hetFs<-which(groups.PhytoMicro[,"group"]=="heterotrophic flagellate")+4 #11 species/taxonomic groups

#ciliates:
ciliates<-which(groups.PhytoMicro[,"group"]=="ciliate")+4 #17 species/taxonomic groups

###EXPLORE DATA & FINAL SELECTION
####Note: the calculated summed biomass carbon are lower than those I calculated for Excel...
#-->It is because I excluded all observations from Ferry Terminal to calculate the summed carbon in R

#Seasonal cycle of the major phytoplankton groups:
#-->2007-2013
#-->diatoms, Phaeocystis sp., Cryptophyceae (only available till 2010), D. speculum, (N. scintillans)
#df<-data.frame(
 # year=datall$year,
 # month=datall$month,
 # date=datall$Date_Time,
 # diatoms=unname(apply(datall[,diatoms],1,sum,na.rm=TRUE)),
 # Phaeocystis_sp=datall[,Phaeo],
 # Cryptophyceae=datall[,Crypto],
 # D.speculum=datall[,Dspec],
 # N.scintillans=datall[,"N._scintillans_C"])
#df<-data.frame(year=rep(df$year,5),date=rep(df$date,5),
              # species=c(rep("diatoms",nrow(df)),rep("Phaeocystis sp.",nrow(df)),
                        # rep("Cryptophyceae",nrow(df)),rep("Dictyocha speculum",nrow(df)),
                        # rep("Noctiluca scintillans",nrow(df))),
              # biomass=c(df$diatoms,df$Phaeocystis_sp,
                        # df$Cryptophyceae,df$D.speculum,df$N.scintillans))
#for(a in 2007:2013){
 # plot<-ggplot(df[which(df[,"year"]==a),],aes(x=date,y=biomass,color=species,shape=species))+
  #  theme_classic()+
  #  geom_point(size=0.7)+
  # geom_smooth(mapping=aes(y=biomass,x=date,color=species),data=df[which(df[,"year"]==a),],method="loess",size=0.5, linetype="solid",se=FALSE)+
  #  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
  #                     breaks=c("diatoms","Phaeocystis sp.","Noctiluca scintillans","Cryptophyceae","Dictyocha speculum"))+
  #  scale_shape_manual("", values=c(15,16,17,18,3),
  #                     breaks=c("diatoms","Phaeocystis sp.","Noctiluca scintillans","Cryptophyceae","Dictyocha speculum"))+
  #  ylab("Biomass carbon [ug C/L]")+
  #  ggtitle(as.character(a))+
  #  theme(plot.title=element_text(face="bold",size=14,color="navy",hjust=0.5),
  #        axis.text.x=element_text(size=6,face="bold",color="black"),
  #        axis.text.y=element_text(face="bold",color="slategray3"),
  #        axis.title.x=element_blank(),
  #        axis.title.y=element_text(face="bold", size=6,color="slategray3"))+
  #  facet_wrap(~species,scales="free",ncol=2)
 # print(plot)}

#############################################################
##DIATOMS########
#Group composition:
carbonSum<-unname(apply(datall[,diatoms],2,FUN=sum, na.rm=TRUE))
sum(carbonSum)#summed carbon biomass of 23252.25 ug C/L
#NA values:
NAs<-unname(apply(datall[,diatoms],2,FUN=function(x)length(which(is.na(x)==TRUE))))
#Potential outtakes:
###OPTION 1####<5% contribution to total diatom biomass
#biom5<-0.05*sum(carbonSum)#1162.612 ug C/L
#carbonSum.rm<-carbonSum
#out<-carbon.excluded<-vector()
#repeat{
#mydiatom<-which(carbonSum.rm==min(carbonSum.rm,na.rm=TRUE))[1]
#if((sum(carbon.excluded)+carbonSum.rm[mydiatom])>biom5){break
#}else{
#    out<-c(out,mydiatom)
#   carbon.excluded<-c(carbon.excluded,carbonSum.rm[mydiatom])
#    carbonSum.rm[mydiatom]<-NA}}
#-->79 species/taxonomic groups form together <5% of total diatom biomass
#-->summed carbon of 1129.062 ug C/L

####OPTION 2###species/taxonomic groups with NA-inflated data sets
#=these groups have continuously not been recorded during 1 year or several years (min of 45 NA values)
#Biomass carbon distribution of not-continuously & continuously recorded spp/tax. groups:
diat0<-which(NAs==0)#32 species/taxonomic groups
diatNA<-which(NAs!=0)#83 species/taxonomic groups
#df<-data.frame(species=c(groups.PhytoMicro[(diatoms-4)[diat0],"species"],groups.PhytoMicro[(diatoms-4)[diatNA],"species"]),
                        #   carbon=c(carbonSum[diat0],carbonSum[diatNA]),NAs=c(NAs[diat0],NAs[diatNA]),
                        #            type=c(rep("continuously recorded",length(diat0)),rep("missing records",length(diatNA))))
#type<-c("continuously recorded","missing records")
#for(i in 1:length(type)){
#plot<-ggplot(df[which(df[,"type"]==type[i]),],aes(x=species,y=carbon))+
#  theme_minimal()+
#  ylab("Summed biomass carbon\n[ug C/L]")+
#  ggtitle(type[i])+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine",color="mediumaquamarine")+
#  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
#        axis.text.x=element_text(angle=90,size=7,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
#  ylim(0,5000)
#print(plot)}
#% contribution to total biomass carbon
#sum(carbonSum[diat0])*100/sum(carbonSum)#spp/tax. groups which have been continuously recorded make up 65.2423% of total biomass carbon
#sum(carbonSum[diatNA])*100/sum(carbonSum)#34.7577%
#########
#Add 3 species with uncontinuous records, 
#but continuously recorded in 2009-2011 & high biomass carbon contribution:
#-->Actinocyclus sp.
diatIN1<-which(colnames(datall[,diatoms])=="Actinocyclus_sp._C")
#-->Mediopyxis helysia
diatIN2<-which(colnames(datall[,diatoms])=="M._helysia_C")
#-->Odontella sinensis
diatIN3<-which(colnames(datall[,diatoms])=="O._sinensis_C")
diatNA<-diatNA[-c(which(diatNA==diatIN1),which(diatNA==diatIN2),
                  which(diatNA==diatIN3))]
out<-diatNA

#####FINAL SELECTION (OPTION 2):35 species/taxonomic groups
diatoms<-diatoms[-out]
sum(carbonSum[-out])#-->summed carbon of 19841.93 ug C/L
df<-data.frame(member=groups.PhytoMicro[(diatoms-4),"species"],
               carbon=carbonSum[-out], NAs=NAs[-out])
#a)Plot the carbon distribution:
#ggplot(df,aes(x=member,y=carbon))+
#  theme_minimal()+
#  ylab("Summed biomass carbon\n[ug C/L]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine",color="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=7,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=8,color="slategray3"))

##############################################################
##Phototrophic flagellates#########
#-->represented by Phaeocystis sp. only
#sum(datall[,Phaeo],na.rm=TRUE)#9701.266 ug C/L

##############################################################
##Dinoflagellates##########
#Group composition:
carbonSum<-unname(apply(datall[,dinos],2,FUN=sum, na.rm=TRUE))
#sum(carbonSum)#summed carbon of 55364.21 ug C/L
#df<-data.frame(member=groups.PhytoMicro[which(groups.PhytoMicro[,"group"]=="dinoflagellate"),"species"],
#               carbon=carbonSum)
#Plot the biomass carbon distribution
#ggplot(df,aes(x=member,y=carbon))+
#  theme_minimal()+
#  ylab("Summed biomass carbon\n[ug C/L]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine",color="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=6,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=8,color="slategray3"))
#NA values:
NAs<-unname(apply(datall[,dinos],2,FUN=function(x)length(which(is.na(x)==TRUE))))
#Potential outtakes:
#<5% contribution to total dinoflagellate biomass
biom5<-0.05*sum(carbonSum)#2768.211 ug C/L
carbonSum.rm<-carbonSum
out<-carbon.excluded<-vector()
repeat{
  mydino<-which(carbonSum.rm==min(carbonSum.rm,na.rm=TRUE))[1]
  if((sum(carbon.excluded)+carbonSum.rm[mydino])>biom5){break
  }else{
    out<-c(out,mydino)
    carbon.excluded<-c(carbon.excluded,carbonSum.rm[mydino])
    carbonSum.rm[mydino]<-NA}}
#-->54 species/taxonomic groups form together <5% of total dinoflagellate biomass
#-->Noctiluca scintillans & Gyrodinium sp.make up together >95% of dinoflagellate biomass
#FINAL SELECTION: only select Noctiluca scintillans 
#-->very different of Gyrodinium (much larger, can be classified as mesozooplankton)
gyrodinium<-which(colnames(datall)=="Gyrodinium_sp._C")
out<-c(out,which(dinos==gyrodinium))
dinos<-dinos[-out]#rm low abundant spp+Gyrodinium
#sum(carbonSum[-out])#-->summed carbon of 52 232.67 ug C/L

######################################################################
##Heterotrophic flagellates############
#Group composition:
#carbonSum<-unname(apply(datall[,hetFs],2,FUN=sum, na.rm=TRUE))
#sum(carbonSum)#summed carbon of 2234.422 ug C/L
#df<-data.frame(member=groups.PhytoMicro[which(groups.PhytoMicro[,"group"]=="heterotrophic flagellate"),"species"],
#               carbon=carbonSum)
#ggplot(df,aes(x=member,y=carbon))+
#  theme_minimal()+
#  ylab("Summed biomass carbon\n[ug C/L]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine",color="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=8,color="slategray3"))
#NA values:
#NAs<-unname(apply(datall[,hetFs],2,FUN=function(x)length(which(is.na(x)==TRUE))))
#Potential outtakes:
#<5% contribution to total het. flagellate biomass
#biom5<-0.05*sum(carbonSum)#111.7211 ug C/L
#carbonSum.rm<-carbonSum
#out<-carbon.excluded<-vector()
#repeat{
#  myhetF<-which(carbonSum.rm==min(carbonSum.rm,na.rm=TRUE))[1]
#  if((sum(carbon.excluded)+carbonSum.rm[myhetF])>biom5){break
#  }else{
#    out<-c(out,myhetF)
#    carbon.excluded<-c(carbon.excluded,carbonSum.rm[myhetF])
#   carbonSum.rm[myhetF]<-NA}}
#-->8 species/taxonomic groups form together <5% of total het. flagellate biomass
#FINAL SELECTION: Ebria tripartita & fractionated/ellipsoid flagellates
#hetFs<-hetFs[-out]
#sum(carbonSum[-out])#summed carbon of 2128.862 ug C/L
#df$NAs<-NAs ; df[-out,c("member","NAs")]#high number of NA values

#########################################################################
##Ciliates############
#Group composition:
carbonSum<-unname(apply(datall[,ciliates],2,FUN=sum, na.rm=TRUE))
#sum(carbonSum)#summed carbon of 2670.617 ug C/L
#df<-data.frame(member=groups.PhytoMicro[which(groups.PhytoMicro[,"group"]=="ciliate"),"species"],
#               carbon=carbonSum)
#ggplot(df,aes(x=member,y=carbon))+
#  theme_minimal()+
#  ylab("Summed biomass carbon\n[ug C/L]")+
#  geom_bar(position = "dodge", stat="identity", width = 0.5, fill="mediumaquamarine", color="mediumaquamarine")+
#  theme(axis.text.x=element_text(angle=90,size=8,face="bold",color="navy"),
#        axis.text.y=element_text(face="bold",color="slategray3"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_text(face="bold", size=8,color="slategray3"))
#NA values:
NAs<-unname(apply(datall[,ciliates],2,FUN=function(x)length(which(is.na(x)==TRUE))))
#Potential outtakes:
####OPTION 1###<5% contribution to total het. flagellate biomass
#biom5<-0.05*sum(carbonSum)#133.5309 ug C/L
#carbonSum.rm<-carbonSum
#out<-carbon.excluded<-vector()
#repeat{
#  mycil<-which(carbonSum.rm==min(carbonSum.rm,na.rm=TRUE))[1]
#  if((sum(carbon.excluded)+carbonSum.rm[mycil])>biom5){break
#  }else{
#    out<-c(out,mycil)
#    carbon.excluded<-c(carbon.excluded,carbonSum.rm[mycil])
#    carbonSum.rm[mycil]<-NA}}
#-->13 species/taxonomic groups form together <5% of total ciliate biomass

###OPTION 2####Only include species/taxonomic groups that were continuously recorded
#-->exclude all groups that were not recorded during 1-2years or more (minimum of 93 NAs)
ciliates0<-which(NAs==0)#4 species/taxonomic groups
ciliatesNA<-which(NAs!=0)#13 species/taxonomic groups
#sum(carbonSum[ciliates0])*100/sum(carbonSum)#selected groups make up 79.68135% of total biomass carbon
out<-ciliatesNA

#FINAL SELECTION (OPTION2): Acineta sp., Laboea strobila & undetermined ciliates/tintinnids
ciliates<-ciliates[-out]
#sum(carbonSum[-out])#summed carbon of 2542.345 ug C/L
#df$NAs<-NAs ; df[-out,c("member","NAs")]#no NA values, except undetermined tintinnids

#####################################################################################
#DATA POOLING
#Pool data of each functional group (=sum up biomass carbon of members):
#-->updated
datall2<-data.frame(
  year=datall$year,
  month=datall$month,
  date=datall$Date_Time,
  diatoms=unname(apply(datall[,diatoms],1,sum,na.rm=TRUE)),
  N.scintillans=datall[,dinos],
  Phaeocystis_sp=datall[,Phaeo],
  #het_flagellates=unname(apply(datall[,hetFs],1,sum,na.rm=TRUE)),
  ciliates=unname(apply(datall[,ciliates],1,sum,na.rm=TRUE)))
group.names<-c("diatoms","N.scintillans", "Phaeocystis_sp", "ciliates")
#write_xlsx(datall2,path="C:/Hannah/Microplankton2007-2013.xlsx")
###################################################################################
###CALCULATE PARAMETER ESTIMATES
source("R functions/get_estimates_MULTIparam.R")
data.PhytoMicro<-get_estimates_MULTIparam(dat=datall2,timeseries=2007:2013,
                                          stations=c("ST1"),colnames.params=group.names)[[2]]
#write_xlsx(data.PhytoMicro,path="C:/Hannah/Microplankton2007-2013_seasonalAVG.xlsx")
# + compare average/median biomass of N. scintillans
source("R functions/get_estimates_ONEparam.AVG_MEDIAN.R")
Noctiluca<-get_estimates_ONEparam.AvgMed(dat=datall2,timeseries=2009:2011,
                        stations=c("ST1"),colname.param="N.scintillans")[[2]]
Noctiluca$season<-paste(Noctiluca$year,Noctiluca$season, sep=" ")
Noctiluca<-Noctiluca[1:which(Noctiluca[,"season"]=="2011 winter"),]
#Export seasonal median Noctiluca biomass as excel file:
#write_xlsx(Noctiluca, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/CBIOM_Nsci.xlsx")

# + compare average/median biomass of diatoms in spring 2009-winter 2010/11:
Diatoms<-get_estimates_ONEparam.AvgMed(dat=datall2,timeseries=2009:2011,
                                    stations=c("ST1"),colname.param="diatoms")[[2]]
Diatoms$season<-paste(Diatoms$year,Diatoms$season, sep=" ")
Diatoms<-Diatoms[2:which(Diatoms[,"season"]=="2011 winter"),]
#Export seasonal median diatom biomass as excel file:
#write_xlsx(Diatoms, 
#path="C:/Hannah/Biological Oceanography/Master Thesis/Project/Model parameters/CBIOM_Dia.xlsx")

###################################################################################
##SEASONAL CYCLE OF THE SELECTED MICROPLANKTON GROUPS (2007-2013)##
response<-paste(group.names,"seasonal_AVG",sep="_")
title<-c("diatoms","N. scintillans","Phaeocystis sp.","ciliates")
for(i in 1:length(response)){
  df<-data.PhytoMicro[,c("year","season",response[i])]
  colnames(df)<-c("year","season","Group")
plot<-ggplot(df,aes(x=year,y=Group,color=season,shape=season))+
  theme_classic()+
  geom_point(size=3)+
  geom_smooth(mapping=aes(y=Group,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon","darkgoldenrod","black"),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  scale_shape_manual("", values=c(15,16,17,18,3),
                     breaks=c("winter","spring","summer","autumn","annual"))+
  ylab("Biomass carbon [ug C/L]")+
  ggtitle(paste("Seasonal cycle of",title[i],"in SRB",sep=" "))+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="slategray3"))
print(plot)}

####################################################################################
setwd(plot.working_dir)
#SEASONALITY OF MICROPLANKTON IN SRB: 04.12.2008-17.02.2011
#-->no sampling in Jan+Feb2010=ice-cold winter
MICRO2009_11<-datall[,c("Event","year","month","Date_Time")]
MICRO2009_11<-cbind(MICRO2009_11,datall[,c(diatoms,dinos,Phaeo,ciliates)])
MICRO2009_11$Diatoms<-apply(datall[,diatoms],1,sum)
MICRO2009_11$Ciliates<-apply(datall[,ciliates],1,sum)
MICRO2009_11<-MICRO2009_11[which(MICRO2009_11[,"Date_Time"]>="2008-12-04"),]
MICRO2009_11<-MICRO2009_11[which(MICRO2009_11[,"Date_Time"]<="2011-02-17"),]

MIC2009_11<-get_estimates_MULTIparam(dat=MICRO2009_11,timeseries=2008:2011,
                                       stations=c("ST1"),
                                       colnames.params=c("Diatoms","Phaeocystis_sp._C","N._scintillans_C","Ciliates"))[[1]]
MIC2009_11<-MIC2009_11[-which(MIC2009_11[,"Nobs_ST1"]==0),]
MIC2009_11$Date<-paste(MIC2009_11$year,MIC2009_11$month,sep="-")
MIC2009_11$Date<-factor(MIC2009_11$Date,levels=MIC2009_11$Date)
#Plot seasonality
mic2009_11<-data.frame(Time=rep(MIC2009_11$Date,4),
                       Comp=c(rep("Diatoms",nrow(MIC2009_11)), rep("Phaeocystis sp.",nrow(MIC2009_11)),
                              rep("N. scintillans",nrow(MIC2009_11)),rep("Ciliates",nrow(MIC2009_11))),
                       Carbon=c(MIC2009_11$Diatoms_monthly_AVG,
                                MIC2009_11$Phaeocystis_sp._C_monthly_AVG,
                                MIC2009_11$N._scintillans_C_monthly_AVG,
                                MIC2009_11$Ciliates_monthly_AVG))
mic2009_11$Time<-factor(mic2009_11$Time,levels=MIC2009_11$Date)
mic2009_11$Comp<-factor(mic2009_11$Comp,levels=c("Diatoms","Phaeocystis sp.","N. scintillans","Ciliates"))
#Store data:
#write.csv(mic2009_11,file=paste(plot.working_dir,"Biomass/microplankton.csv",sep=""))
#open graphical device:
#-->1.5 column width
pdf(file="Microplankton.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(mic2009_11,aes(x=Time,y=Carbon,group=Comp,color=Comp))+
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

######DIATOMS:
dia.names<-c()
for(i in 5:39){
dia.name<-colnames(MICRO2009_11)[i]
if(dia.name=="O._rhombus_f_trigona_C")dia.name<-"O._rhombus.ft_C"
mydia<-paste(strsplit(dia.name,split="_",fixed=TRUE)[[1]][1],
                 strsplit(dia.name,split="_",fixed=TRUE)[[1]][2],sep=" ")
dia.names<-c(dia.names,mydia)}

DIA2009_11<-MICRO2009_11[,c(1:39)]
colnames(DIA2009_11)[5:39]<-dia.names
DIA2009_11<-get_estimates_MULTIparam(dat=DIA2009_11,timeseries=2008:2011,
                           stations=c("ST1"),
                           colnames.params=dia.names)[[1]]
DIA2009_11<-DIA2009_11[-which(DIA2009_11[,"Nobs_ST1"]==0),]
DIA2009_11$Date<-paste(DIA2009_11$year,DIA2009_11$month,sep="-")
DIA2009_11$Date<-factor(DIA2009_11$Date,levels=DIA2009_11$Date)
dia2009_11<-as.data.frame(matrix(NA,ncol=3,nrow=0,dimnames=list(NULL,c("Time","Spp","Carbon"))))
for(i in 1:length(dia.names)){
  col.mydia<-paste(dia.names[i],"monthly_AVG",sep="_")
  mydia<-data.frame(Time=DIA2009_11$Date,Spp=rep(dia.names[i],nrow(DIA2009_11)),
                              Carbon=DIA2009_11[,col.mydia])
  dia2009_11<-rbind(dia2009_11,mydia)}
dia2009_11$Time<-factor(dia2009_11$Time,levels=DIA2009_11$Date)
cols<-rep("powderblue",length(unique(dia2009_11$Spp)))
cols[which(unique(dia2009_11$Spp)=="O. aurita")]<-"hotpink"
cols[which(unique(dia2009_11$Spp)=="B. brockmannii")]<-"mediumpurple"
cols[which(unique(dia2009_11$Spp)=="R. imbricata")]<-"lightseagreen"
cols[which(unique(dia2009_11$Spp)=="Chaetoceros sp.")]<-"violetred"
cols[which(unique(dia2009_11$Spp)=="P-n delicatissima")]<-"olivedrab"
cols[which(unique(dia2009_11$Spp)=="L. danicus")]<-"cyan"
cols[which(unique(dia2009_11$Spp)=="O. mobiliensis")]<-"orange"
cols[which(unique(dia2009_11$Spp)=="O. rhombus")]<-"black"
cols[which(unique(dia2009_11$Spp)=="A. glacialis")]<-"red"
cols[which(unique(dia2009_11$Spp)=="Actinocyclus sp.")]<-"royalblue"
cols[which(unique(dia2009_11$Spp)=="M. helysia")]<-"black"
cols[which(unique(dia2009_11$Spp)=="O. sinensis")]<-"yellow"

#open graphical device:
#-->1.5 column width
pdf(file="Diatoms.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(dia2009_11,aes(x=Time,y=Carbon,group=Spp,color=Spp))+
  theme_minimal()+
  geom_point()+
  geom_line()+
  ylab("mg C/m^3")+
  scale_color_manual("",values=cols,breaks=unique(dia2009_11$Spp))+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="black"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="black"))
# close the graphical device:
dev.off() 

##################################################################
#MEAN CARBON CONTRIBUTION (based on average of contribution of each species/group to respective carbon pool at each sampling day)
#meanCC=mean(carbon of species i at sampling day x/carbon pool at sampling day x)
#-->See method section about "Weighted Community Biomass" of Julien Meunier
##################################################################
setwd(working_dir)
source("R functions/carbon_contributions.R")

#A)####CARBON CONTRIBUTION OF EACH INCLUDED SPP TO COMPARTMENT CARBON
##############ACROSS SEASONS 2009/2010 (including Dec 2008 & Jan.+Feb. 2011)
micro.species<-micro.cols<-list()
micro.comps<-c("diatoms","phytoplankton","ciliates")
#-->diatom species to diatoms
micro.species[[1]]<-dia.names
micro.cols[[1]]<-colnames(datall)[diatoms]
#-->diatoms to total phytoplankton
micro.species[[2]]<-c("diatoms","Phaeocystis sp.")
micro.cols[[2]]<-c("Diatoms","Phaeocystis_sp._C")
#-->ciliate species to ciliates
micro.species[[3]]<-spp.PhytoMicro[(ciliates-4),"species"]
micro.cols[[3]]<-colnames(datall)[ciliates]

#Extract df with total microzooplankton carbon used in the model
TOTCARB.m<-MICRO2009_11[which(MICRO2009_11[,"Date_Time"]>="2008-12-04"),]
TOTCARB.m<-TOTCARB.m[which(TOTCARB.m[,"Date_Time"]<="2011-02-17"),]
colnames(TOTCARB.m)[4]<-"date"

CARBcontr.MICRO<-spptocompC(dat=TOTCARB.m,stations="1",species=micro.species,
                            columns=micro.cols,comps=micro.comps)

#B)####CARBON CONTRIBUTION OF EACH COMPARTMENT TO TOTAL FUNC. GROUP CARBON
##############ACROSS SEASONS 2009/2010 (including Dec 2008 & Jan.+Feb. 2011)
#-->diatoms to total diatoms & total phytoplankton
#-->Phaeocystis to total phytoplankton
#-->phytoplankton to total phytoplankton
#-->ciliates to total ciliates & total protozooplankton (het. flagellates, ciliates, dinoflagellates, except N.scintillans)
#-->dinoflagellates (not represented in model) to total protozooplankton
#-->het. flagellates (not represented in model) to total protozooplankton
#-->model microplankton to total microplankton (except Noctiluca scintillans)

#Extract df with total microplankton carbon of SRB: Dec. 2008-Feb.2011
TOTCARB<-datall[which(datall[,"Date_Time"]>="2008-12-04"),]
TOTCARB<-TOTCARB[which(TOTCARB[,"Date_Time"]<="2011-02-17"),]
colnames(TOTCARB)[3]<-"date"

#Import df with grouping variables of all SRB microplankton species
SRB.spp<-fread(file = "Phyto-_Microplankton.spp_J.Rick/SRB_Total_Microplankton_Carbon.csv", na.strings = "", dec = "," , data.table = FALSE)
#Reorder df according to TOTCARB column order
order<-data.frame(id=1:224,column=colnames(TOTCARB)[5:228])
SRB.spp<-merge(SRB.spp,order,by="column")
SRB.spp<-SRB.spp[order(SRB.spp$id,decreasing=FALSE),]

#Prepare df containing model + total SRB carbon:

#Total model phytoplankton carbon:
TOTCARB.m$Phytoplankton<-apply(TOTCARB.m[,c("Diatoms","Phaeocystis_sp._C")],1,sum,na.rm=TRUE)
#Total model microplankton carbon:
TOTCARB.m$Microplankton<-apply(TOTCARB.m[,c("Diatoms","Phaeocystis_sp._C","Ciliates")],1,sum,na.rm=TRUE)

#Total SRB diatom carbon:
TOTCARB$Diatoms.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"group"]=="diatom"))],1,sum,na.rm=TRUE)
#Total SRB ciliate carbon:
TOTCARB$Ciliates.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"group"]=="ciliate"))],1,sum,na.rm=TRUE)
#Total SRB hetF carbon:
TOTCARB$HetFs.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"group"]=="heterotrophic flagellate"))],1,sum,na.rm=TRUE)
#Total SRB dinoflagellate carbon:
TOTCARB$Dinos.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"group"]=="dinoflagellate"))],1,sum,na.rm=TRUE)
#Total SRB phytoplankton carbon:
TOTCARB$Phytoplankton.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"funcGroup"]=="phytoplankton"))],1,sum,na.rm=TRUE)
#Total SRB protozooplankton carbon:
TOTCARB$Protozooplankton.SRB<-apply(TOTCARB[,(4+which(SRB.spp[,"funcGroup"]=="protozooplankton"))],1,sum,na.rm=TRUE)
#Total SRB microzooplankton carbon:
TOTCARB$Microplankton.SRB<-apply(TOTCARB[,c("Phytoplankton.SRB","Protozooplankton.SRB")],1,sum,na.rm=TRUE)

TOTCARB.SRB<-merge(TOTCARB.m[,c("Event","year","month","date","Diatoms","Phaeocystis_sp._C","Ciliates",
                          "Phytoplankton","Microplankton")],
             TOTCARB[,c("Event","year","month","date","Diatoms.SRB","Ciliates.SRB","HetFs.SRB","Dinos.SRB",
                        "Phytoplankton.SRB","Protozooplankton.SRB","Microplankton.SRB")], by=c("Event","year","month","date"))

comps<-c(rep("diatoms",2),"Phaeocystis sp.","diatoms + Phaeocystis sp.",
         rep("ciliates",2),"dinoflagellates","heterotrophic flagellates",
         "diatoms + Phaeocystis sp. + ciliates")
cols.comps<-c(rep("Diatoms",2),"Phaeocystis_sp._C","Phytoplankton",
              rep("Ciliates",2),"Dinos.SRB","HetFs.SRB",
              "Microplankton")
totC<-c("SRB diatoms",rep("SRB phytoplankton",3),
         "SRB ciliates",rep("SRB protozooplankton",3),"SRB microplankton")
cols.totC<-c("Diatoms.SRB",rep("Phytoplankton.SRB",3),
        "Ciliates.SRB",rep("Protozooplankton.SRB",3),"Microplankton.SRB")

CARBcontr.MICRO_SRB<-comptototC(dat=TOTCARB.SRB,stations="1",comps,
                            cols.comps,totC,cols.totC)

#Save results:
setwd(plot.working_dir)
write_xlsx(CARBcontr.MICRO,"Microplankton.xlsx")
write_xlsx(CARBcontr.MICRO_SRB,"Microplankton_SRB.xlsx")

#Plot results:
#open graphical device:
#-->1.5 column width
pdf(file="MicroplanktonSpeciesComposition.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
plot(plot.Ccontributions(CARBcontr.MICRO))
# close the graphical device:
dev.off() 

#open graphical device:
#-->1.5 column width
pdf(file="TotalSRBMicroplankton.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
plot(plot.Ccontributions(CARBcontr.MICRO_SRB[-which(CARBcontr.MICRO_SRB[,"Compartment"]=="diatoms + Phaeocystis sp."),]))
#-->summed protozooplankton members often <100% because some 
#   protozooplankton were not classified in any of the membergroups (eg. cyanobacteria)
# close the graphical device:
dev.off()
