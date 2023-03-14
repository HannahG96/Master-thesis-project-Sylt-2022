##### ADDITIONAL PLOTS #####

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/")

#load packages:
library("data.table")
library("ggplot2")
library("writexl")
library("LIM")

#Compartments:
comps<-c("Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Bac","Doc","Poc")
n<-length(comps)
living<-c("Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Bac")
nonliving<-c("Doc","Poc")
funcGroups<-c("Diatoms","Phaeocystis sp.","Ciliates","Noctiluca sp.",rep("Mesozooplankton",3),
              rep("Meroplankton",3),rep("Jellyfish",4),"Herring","Bacteria",rep("Detritus",2))
liv<-16
nl<-2
externals<-c("Imp","Resp","Exp")
seasons<-c("2009 spring","2009 summer","2009 autumn","2010 spring","2010 summer","2010 autumn")
#
#
#
#
##### SEASONAL BIOMASS OF FUNCTIONAL GROUPS (LOG-SCALE) #####
#Import information on living compartment biomass:
CBIOM<-fread(file = "Parameters/CBIOM.csv", na.strings = "", dec = "," , data.table = FALSE)
CBIOM<-CBIOM[-c(which(CBIOM$season=="2010 winter"),which(CBIOM$season=="2011 winter")),]
CBIOM$Bac<-CBIOM$Doc<-CBIOM$Poc<-NA
CBIOM<-CBIOM[,c("year","season",comps)] #reorder biomass data according to compartment order
BIOM<-as.matrix(CBIOM[,3:ncol(CBIOM)])
colnames(BIOM)<-funcGroups
rownames(BIOM)<-seasons
BIOM[which(BIOM==0)]<-0.001 #attribute absent compartments a minimum biomass value
#
#
groups<-c("Diatoms","Phaeocystis sp.","Ciliates","Noctiluca sp.","Mesozooplankton","Meroplankton",
          "Jellyfish","Herring")
funcBIOM<-as.data.frame(matrix(0,nrow=0,ncol=3,dimnames=list(NULL,c("season","Group","Biomass"))))
for(i in 1:length(groups)){
  if(length(which(funcGroups==groups[i]))>1){
  seasonal.biomass<-as.numeric(apply(as.data.frame(BIOM[,which(colnames(BIOM)==groups[i])]),1,sum))
  }else{seasonal.biomass<-as.numeric(BIOM[,groups[i]])}
  dat<-data.frame(season=seasons,Group=rep(groups[i],length(seasons)),Biomass=seasonal.biomass)
  funcBIOM<-rbind(funcBIOM,dat)}
funcBIOM$season<-factor(funcBIOM$season,levels=seasons)
funcBIOM$Group<-factor(funcBIOM$Group,levels=groups)
funcBIOM$YEAR<-funcBIOM$SEASON<-NA
for(i in 1:nrow(funcBIOM)){
  funcBIOM[i,"YEAR"]<-strsplit(as.character(funcBIOM[i,"season"]),split=" ",fixed=TRUE)[[1]][1]
  funcBIOM[i,"SEASON"]<-strsplit(as.character(funcBIOM[i,"season"]),split=" ",fixed=TRUE)[[1]][2]}
YEARs<-c(2009,2010)
SEASONs<-c("2009 spring","2010 spring","2009 summer","2010 summer","2009 autumn","2010 autumn")
funcBIOM$season<-factor(funcBIOM$season,levels=SEASONs)
funcBIOM$SEASON<-factor(funcBIOM$SEASON,levels=c("spring","summer","autumn"))
funcBIOM$YEAR<-factor(funcBIOM$YEAR, levels=YEARs)

#Plot data:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/")
#open graphical device:
#-->1.5 column width
pdf(file="Figures/Biomass.pdf",         # File name
   width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(funcBIOM, aes(x = season, y = Biomass)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(color=SEASON,shape=YEAR)) +
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(16,15),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("Biomass [mg C/m3]")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~Group,ncol=2,scales="free")
# close the graphical device:
dev.off() 

## Alternative plots ##
#-->seasonal dynamics of each functional group
plot.working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Community composition/Biomass/"
setwd(plot.working_dir)

#Import biomass data:
files<-paste(c("microplankton","zooplankton","meroplankton","jellyfish"),".csv",sep="")
Biomass<-as.data.frame(matrix(NA,nrow=0,ncol=3,dimnames=list(NULL,c("Time","Comp","Carbon"))))
for(i in 1:length(files)){
  dat<-fread(file = files[i], na.strings = "NA", dec = "," , data.table = FALSE,header=TRUE)[2:4]
  Biomass<-rbind(Biomass,dat)}
#Re-organise data:
#1)year, month & season column
Biomass$Year<-Biomass$Month<-NA
for(i in 1:nrow(Biomass)){
  Biomass[i,"Year"]<-as.numeric(strsplit(Biomass[i,"Time"],split="-")[[1]][1])
  Biomass[i,"Month"]<-as.numeric(strsplit(Biomass[i,"Time"],split="-")[[1]][2])}
monthLabs<-c("Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec.")
seasonLabs<-c(rep("Winter",2),rep("Spring",3),rep("Summer",3),rep("Autumn",3),"Winter")
Biomass$monthLabs<-Biomass$seasonLabs<-NA
for(i in 1:12){
  Biomass[which(Biomass$Month==i),"monthLabs"]<-monthLabs[i]
  Biomass[which(Biomass$Month==i),"seasonLabs"]<-seasonLabs[i]}
#2)compartment labels
compLabs<-c("Dia","Phae","Cil","Nsci","Cop","Clado","Tun","Biv","Gastr","Poly","Hydro","Mlei","Ppil","Bcu")
compNames<-c("Diatoms",         "Phaeocystis sp.",   "Ciliates",        "N. scintillans",     
            "Copepods",        "Cladocerans",  "tunicates" ,     "Bivalvia",  "gastropods",  "Polychaetes",    
            "Hydromedusae",  "M.leidyi",        "P.pileus",        "B.cucumis" )
Biomass$compLabs<-NA
for(i in 1:length(compLabs)){
  Biomass[which(Biomass$Comp==compNames[i]),"compLabs"]<-compLabs[i]}
#3)functional groups
compGroups<-c(rep("Phytoplankton",2),"Ciliates","Noctiluca sp.",rep("Mesozooplankton",3),
              rep("Meroplankton",3),rep("Jellyfish",4))
Biomass$compGroups<-NA
for(i in 1:length(compGroups)){
  Biomass[which(Biomass$compLabs==compLabs[i]),"compGroups"]<-compGroups[i]}
#4)re-class columns
Biomass$monthLabs<-factor(Biomass$monthLabs,levels=monthLabs)
Biomass$compLabs<-factor(Biomass$compLabs,levels=compLabs)
Biomass$compGroups<-factor(Biomass$compGroups,levels=unique(compGroups))
Biomass$Carbon<-as.numeric(Biomass$Carbon)
#4)remove years 2008 & 2011, and winter season
Bio<-Biomass[-c(which(Biomass$Year==2008),which(Biomass$Year==2011)),]
Bio<-Bio[-which(Bio$seasonLabs=="Winter"),]
Bio$Year<-factor(Bio$Year,levels=c(2009,2010))
Bio$seasonLabs<-factor(Bio$seasonLabs,levels=c("Spring","Summer","Autumn"))
#5)classify data
Bio.micro<-Bio[c(which(Bio$compGroups=="Phytoplankton"),which(Bio$compGroups=="Ciliates"),
                 which(Bio$compGroups=="Noctiluca sp.")),]
Bio.zoo<-Bio[c(which(Bio$compGroups=="Mesozooplankton"),which(Bio$compGroups=="Meroplankton")),]
Bio.jelly<-Bio[which(Bio$compGroups=="Jellyfish"),]

### Plot data:
#microplankton
pdf(file="Biomass_micro.pdf",         # File name
    width = 6, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(Bio.micro,aes(x=monthLabs,y=Carbon))+
  theme_bw()+
  geom_line(aes(group=Year,linetype=Year),size=0.5)+ 
  geom_point(size=2,aes(shape=Year,color=seasonLabs)) +
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("Spring","Summer","Autumn"))+
  xlab("")+
  ylab(expression("mg C/m" ^ "3"))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_text(face="bold",size=7,angle=35),
        axis.title.y = element_text(face="bold", size=10, color="black"),
        #aspect.ratio = 0.35,
        strip.background = element_blank(),
        strip.text=element_text(face="bold"))+
  facet_wrap(~compLabs,ncol=2,scales="free")
# close the graphical device:
dev.off()

#zooplankton
pdf(file="Biomass_zoo.pdf",         # File name
    width = 7.5, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(Bio.zoo,aes(x=monthLabs,y=Carbon))+
  theme_bw()+
  geom_line(aes(group=Year,linetype=Year),size=0.5)+ 
  geom_point(size=2,aes(shape=Year,color=seasonLabs)) +
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("Spring","Summer","Autumn"))+
  xlab("")+
  ylab(expression("mg C/m" ^ "3"))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_text(face="bold",size=7,angle=35),
        axis.title.y = element_text(face="bold", size=10, color="black"),
        #aspect.ratio = 0.35,
        strip.background = element_blank(),
        strip.text=element_text(face="bold"))+
  facet_wrap(~compLabs,ncol=3,scales="free")
# close the graphical device:
dev.off()

#jellyfish
pdf(file="Biomass_jelly.pdf",         # File name
    width = 6, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(Bio.jelly,aes(x=monthLabs,y=Carbon))+
  theme_bw()+
  geom_line(aes(group=Year,linetype=Year),size=0.5)+ 
  geom_point(size=2,aes(shape=Year,color=seasonLabs)) +
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("Spring","Summer","Autumn"))+
  xlab("")+
  ylab(expression("mg C/m" ^ "3"))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_text(face="bold",size=7,angle=35),
        axis.title.y = element_text(face="bold", size=10, color="black"),
        #aspect.ratio = 0.35,
        strip.background = element_blank(),
        strip.text=element_text(face="bold"))+
  facet_wrap(~compLabs,ncol=2,scales="free")
# close the graphical device:
dev.off()

#all groups
pdf(file="Biomass_ALL.pdf",         # File name
    width = 9, height = 6, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(Bio,aes(x=monthLabs,y=Carbon))+
  theme_bw()+
  geom_line(aes(group=Year,linetype=Year),size=0.5)+ 
  geom_point(size=2,aes(shape=Year,color=seasonLabs)) +
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("Spring","Summer","Autumn"))+
  xlab("")+
  ylab(expression("mg C/m" ^ "3"))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_text(face="bold",size=7,angle=35),
        axis.title.y = element_text(face="bold", size=10, color="black"),
        #aspect.ratio = 0.35,
        strip.background = element_blank(),
        strip.text=element_text(face="bold"))+
  facet_wrap(~compLabs,ncol=4,scales="free")
# close the graphical device:
dev.off()
