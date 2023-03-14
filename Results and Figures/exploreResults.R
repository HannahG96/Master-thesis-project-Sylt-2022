######## EXPLORE ENA RESULTS ############

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/")

#load packages:
library("LIM")
library("data.table")
library("limSolve")
library("writexl")
library("ggplot2")

#Compartments:
comps<-c("Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Bac","Doc","Poc")
n<-length(comps)
living<-c("Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Bac")
nonliving<-c("Doc","Poc")
liv<-16
nl<-2
externals<-c("Imp","Resp","Exp")

#Import ENA results and merge them according to season:
seasons<-c("2009 spring","2009 summer","2009 autumn",
           "2010 spring","2010 summer","2010 autumn")
seasonfolders<-c("Spring 2009","Summer 2009","Autumn 2009","Spring 2010","Summer 2010","Autumn 2010")
files<-paste(c("INDS_CI95","ATTR_CI95"),".csv",sep="")
indices<-c("INTER","IMP","EXP","TSTp","TSTf","INTER/TST","IMP/TST","EXP/TST","NPP","NPP/RESP",
           "Herb","Det","Bact","APL","D/H","TE12","TE23","MTL","SOI",
           "H","Ov/DC","A/DC","R/DC","AMI","Ai/DCi","ELD","Ceff","LDq",
           "FCI","NCycles","C_cycled","C_cycled.small","C_cycled.big","GPP","BP","Egestion","Respiration")
INDS_allseasons<-as.data.frame(matrix(NA,nrow=0,ncol=7,dimnames=list(NULL,
                                  c("season","Index","Mean","Q0","Q2.5","Q97.5","Q100"))))
ALL<-list(INDS_allseasons)
for(i in 1:length(seasons)){
  for(u in 1:1){#length(files)
    dat<-fread(file = paste("Results",seasonfolders[i],files[1],sep="/"), na.strings = "", dec = "," , 
               data.table = FALSE)
    datt<-fread(file = paste("Results",seasonfolders[i],files[2],sep="/"), na.strings = "", dec = "," , 
                data.table = FALSE)
    dat<-rbind(dat,datt)
    dat<-cbind(data.frame(season=rep(seasons[i],nrow(dat))),dat)
    ALL[[u]]<-rbind(ALL[[u]],dat)}}
INDS_allseasons<-ALL[[1]]

#### PLOT RESULTS ####
######## COMPARISONS ACROSS SEASONS AND YEARS ########
INDS_allseasons$YEAR<-INDS_allseasons$SEASON<-NA
for(i in 1:nrow(INDS_allseasons)){
  INDS_allseasons[i,"YEAR"]<-strsplit(INDS_allseasons[i,"season"],split=" ",fixed=TRUE)[[1]][1]
  INDS_allseasons[i,"SEASON"]<-strsplit(INDS_allseasons[i,"season"],split=" ",fixed=TRUE)[[1]][2]}
YEARs<-c(2009,2010)
SEASONs<-c("2009 spring","2010 spring","2009 summer","2010 summer","2009 autumn","2010 autumn")

select.inds1<-c("IMP","TSTp","EXP")#Size & activity of the system
select.inds1b<-c("TE12")#Transfer efficiency
select.inds2<-c("H","ELD")#Flow Diversity & Connectivity
select.inds3<-c("A/DC","Ov/DC","R/DC","Ai/DCi")#Flow organisation
select.inds4<-c("Ai/DCi","AMI")#Inherent flow organisation
select.inds5<-c("APL","D/H","FCI")#Retention of energy & recycling

### Size & activity of the system
labs<-c("Total Import","Total System Throughput","Total Export")
for(i in 1:length(select.inds1)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds1[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds1[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds1, labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/SizeActivity.pdf",         # File name
    width = 9, height = 2.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean,ymin=Q2.5, ymax=Q97.5)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                    breaks=c("spring","summer","autumn"),
                  labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                     breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab(expression(paste("mg C/m" ^ "3","/day",sep="")))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="black"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

labs<-c("Transfer Efficiency: TL1->TL2")
for(i in 1:length(select.inds1b)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds1b[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds1b[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds1b, labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/TransferEfficiency.pdf",         # File name
    width = 5.5, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean*100,ymin=Q2.5*100, ymax=Q97.5*100)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("%")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

#Flow diversity & connectivity
labs<-c("Flow Diversity","Effective Link Density")
for(i in 1:length(select.inds2)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds2[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds2[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds2,labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/FlowDiversityConnectivity.pdf",         # File name
    width = 8, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean,ymin=Q2.5, ymax=Q97.5)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("[no unit]")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="black"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

#Flow Organisation
labs<-c("Ascendency","Overhead","Redundancy","Internal Ascendency")
for(i in 1:length(select.inds3)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds3[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds3[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds3,labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/FlowOrganisation.pdf",         # File name
    width = 7, height = 4.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean*100,ymin=Q2.5*100, ymax=Q97.5*100)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("%")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="black"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

#Inherent flow organisation
labs<-c("Internal Ascendency","Average Mutual Information")
for(i in 1:length(select.inds4)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds4[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds4[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds4,labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/InternalFlowOrganisation.pdf",         # File name
    width = 8, height = 3.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean,ymin=Q2.5, ymax=Q97.5)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("[no unit]")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

#Retention of energy & Recycling
labs<-c("Average Path Length","Detritivory-Herbivory Ratio","Finn's Cycling Index")
for(i in 1:length(select.inds5)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds5[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.inds5[i]),])}}
plot.data$Index<-factor(plot.data$Index,levels=select.inds5,labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="Figures/RetentionRecycling.pdf",         # File name
    width = 9, height = 3, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean,ymin=Q2.5, ymax=Q97.5)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("[no unit]")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 

### FINAL PLOT(s) + TABLE(s) ###

#1) Table of indices
ind<-c("IMP","TSTp","EXP","GPP","C_cycled","APL","TE12","D/H","SOI","FCI","H","ELD","A/DC","Ov/DC","R/DC","Ai/DCi","AMI")
Inds<-data.frame(Index=ind)
for(i in 1:length(seasons)){
  myinds<-INDS_allseasons[which(INDS_allseasons$season==seasons[i],),]
  ii<-as.data.frame(matrix(NA,nrow=0,ncol=3,dimnames=list(NULL,c("Index",paste(c("Min","Max"),i,sep="")))))
  for(u in 1:length(ind)){
    myii<-round(myinds[which(myinds$Index==ind[u]),c("Q2.5","Q97.5")],digits=3)
    ii<-rbind(ii,myii)}
 Inds<-cbind(Inds,ii) }
write_xlsx(Inds,path="C:/Hannah/Biological Oceanography/Master Thesis/Project/RESULTS/Indices.xlsx")

#2) Figure of significant indice shifts
#select.indsFINAL<-c("APL","TE12","D/H",
 #                   "ELD","H","SOI",
  #                  "A/DC","Ov/DC","R/DC",
  #                  "Ai/DCi","AMI","FCI")
#labs<-c("Average Path Length","Transfer Efficiency 1->2","Detritivory-Herbivory Ratio",
      #  "Effective Link Density","Flow Diversity","System Omnivory Index",
      # "Ascendency", "Overhead","Redundancy",
      #  "Internal Ascendency","Average Mutual Information","Finn's Cycling Index")
select.indsFINAL<-c("APL","TE12","H",
                    "A/DC","Ov/DC","R/DC",
                    "Ai/DCi")
labs<-c("Average Path Length","Transfer Efficiency 1->2","Flow Diversity",
        "Ascendency", "Overhead","Redundancy",
        "Internal Ascendency")
pct<-c(2,4:7)
for(i in 1:length(select.indsFINAL)){
  if(i==1){
    plot.data<-INDS_allseasons[which(INDS_allseasons[,"Index"]==select.indsFINAL[i]),]
  }else{
    plot.data<-rbind(plot.data,INDS_allseasons[which(INDS_allseasons[,"Index"]==select.indsFINAL[i]),])}}
for(i in pct){
  plot.data[which(plot.data$Index==select.indsFINAL[i]),3:7]<-plot.data[which(plot.data$Index==select.indsFINAL[i]),3:7]*100}
plot.data$Index<-factor(plot.data$Index,levels=select.indsFINAL,labels=labs)
plot.data$season<-factor(plot.data$season,levels=SEASONs)
plot.data$SEASON<-factor(plot.data$SEASON,levels=c("spring","summer","autumn"))
plot.data$YEAR<-factor(plot.data$YEAR,levels=YEARs)
#open graphical device:
#-->1.5 column width
pdf(file="C:/Hannah/Biological Oceanography/Master Thesis/Project/RESULTS/Figures/Indices.pdf",         # File name
    width = 8, height = 6, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(plot.data, aes(x = season, y = Mean,ymin=Q2.5, ymax=Q97.5)) + 
  theme_bw()+
  geom_line(aes(group=YEAR,linetype=YEAR))+ #aes(group=1)
  geom_point(size=2,aes(shape=YEAR)) +
  geom_errorbar(aes(color=SEASON),position=position_dodge(.9), width=0.2)+
  scale_color_manual("Season",values=c("darkcyan","yellowgreen","darkgoldenrod"),
                     breaks=c("spring","summer","autumn"),
                     labels=c("Spring","Summer","Autumn"))+
  scale_linetype_manual("Year",values=c("solid","dotted"),
                        breaks=c(2009,2010))+
  scale_shape_manual("Year",values=c(1,0),
                     breaks=c(2009,2010))+
  xlab("")+
  ylab("")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
  facet_wrap(~Index,scales="free")
# close the graphical device:
dev.off() 
