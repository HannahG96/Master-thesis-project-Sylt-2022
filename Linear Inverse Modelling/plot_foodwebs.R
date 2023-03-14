#### PLOT SEASONAL FOOD WEB MODELS ####
#->based on central value of flows

#set working directory
limwd<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/"
resultwd<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Results/"
setwd(limwd)

#load packages:
library("DiagrammeR")
library("ggplot2")
library("LIM")
library("data.table")

#get functions for plotting:
source("ENA PACKAGE/cycling_analysis.R")
source("ENA PACKAGE/adapted_ENAfuncs.R")
source("ENA PACKAGE/plot_network.R")

columns<-c("2009 spring","2009 summer","2009 autumn","2010 spring","2010 summer","2010 autumn")
seasons<-c("SPRING2009","SUMMER2009","AUTUMN2009","SPRING2010","SUMMER2010","AUTUMN2010")
folders<-c("Spring 2009", "Summer 2009", "Autumn 2009", "Spring 2010", "Summer 2010", "Autumn 2010")
folders.lim<-paste(folders,"Version 5/",sep="/")
outputnames<-paste(seasons,".xlsx",sep="")

#Compartments:
comps<-c("Bac","Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Doc","Poc")
externals<-c("Imp","Resp","Exp")
nl<-2
liv<-16
funcGroups<-c(rep("bottom",3),"mic",rep("meso",4),rep("mero",3),rep("top",5),rep("nl",2))

#####################################################
#Get mean values of flows for each seasonal model:
#####################################################
FW.plots<-list()

networks<-list()
n<-length(comps)
for(i in 1:length(seasons)){
  setwd(limwd)
  readLIM<-Read(paste(folders.lim[i],seasons[i],"_version5.input",sep=""))
  LIM<-Setup(readLIM) 
  ### Create vectors to identify incoming/outgoing compartments:
  FROM<-c()
  TO<-c()
  for(u in 1:length(LIM$Unknowns)){
    from<-strsplit(LIM$Unknowns[u],split="->",fixed=TRUE)[[1]][1]
    FROM<-c(FROM,from)
    to<-strsplit(LIM$Unknowns[u],split="->",fixed=TRUE)[[1]][2]
    TO<-c(TO,to) }  
  ### Get mean value of each flow:
  setwd(resultwd)
  FlowRanges<-fread(file = paste(folders[i],"FLOWS_CI95.csv",sep="/"), na.strings = "", dec = "," , data.table = FALSE)
  ### Produce Transfer matrix and import/export vectors:
  Tstar<-matrix(0,nrow=(length(comps)+3),ncol=(length(comps)+3),
                dimnames=list(c(comps,externals),c(comps,externals)))
  for(u in 1:length(LIM$Unknowns)){
    myrow<-which(rownames(Tstar)==FROM[u])
    mycol<-which(colnames(Tstar)==TO[u])
    Tstar[myrow,mycol]<-FlowRanges[u,"Mean"]}
  Tstar<-Tstar[-c((n+2),(n+3)),-c(n+1)]
  networks[[i]]<-list(Tstar[1:n,1:n],Tstar[n+1,1:n],Tstar[1:n,n+1],Tstar[1:n,n+2])
  names(networks[[i]])<-c("T matrix","Z","R","E")}
names(networks)<-seasons

### Get seasonal TPs of compartments:
TPs_allseasons<-matrix(NA,ncol=length(comps),nrow=length(seasons),dimnames=list(seasons,comps))
for(i in 1:length(networks)){
  TPs_allseasons[i,]<-TP.CTA(inp=networks[[i]][[2]],inter=networks[[i]][[1]],outt=networks[[i]][[4]],
                             diss=networks[[i]][[3]],nl=nl)[[2]]}
max(TPs_allseasons)
breaksY<-c(1,2,2.01,2.1,2.2,2.25,3,3.25)
valsY<-c(0,3,4,5,6,7,8,8.5)

### Get seasonal biomass of compartments
#Import information on living compartment biomass:
setwd(limwd)
CBIOM<-fread(file = "Parameters/CBIOM.csv", na.strings = "", dec = "," , data.table = FALSE)
breaksBIOM<-c(0,0.5,1,2.5,5,10,25,50,75,100,150,200,250)
valsBIOM<-c(0.1,1,2,3,4,4.5,5,6,7,8,9,10,12)

### Define x-axis position of nods:
coordX<-c(12,1.5,4.5,9,6,4.5,0.75,2.25,0,1.5,3,7.5,9,10.5,6,12,7.5,10.5)

### Define additional graph attributes:
attr<-data.frame(attr=c("penwidth","fontsize"),value=c(3,10),attr_type=rep("node",2))

#Plot networks:
for(z in 1:length(seasons)){
  myseason<-columns[z]
  
  ####Define nod size scale based on biomass:
  BIOM<-CBIOM[which(CBIOM[,"season"]==myseason),3:ncol(CBIOM)]
  BIOM$Bac<-BIOM$Doc<-mean(as.numeric(BIOM[1,]))
  #reorder biomass vector according to compartment order
  BIOM<-BIOM[,comps]
  ABSENTcomps<-comps[which(BIOM[1,]==0)]#check absence of compartments 
  BIOM[,ABSENTcomps]<-0.001 #attribute absent compartments a minimum biomass value
  #make biomass vector class numeric:
  BIOM<-as.numeric(BIOM)
  #define biomass scale
  coefBIOM<-rep(0,length(comps))
  for(k in 1:length(BIOM))coefBIOM[which(BIOM>=breaksBIOM[k])]<-valsBIOM[k]
  
  ####Define Y-axis scale based TPs, except non-living+primary-producing compartments:
  coefY<-rep(0,length(comps))
  for(k in 1:length(coefY))coefY[which(TPs_allseasons[z,]>=breaksY[k])]<-valsY[k]
  
  ####Define color code of nods:
  nodcols<-c("skyblue3",rep("mediumseagreen",2),"mediumvioletred",rep("tomato1",4),rep("yellow3",3),
             rep("turquoise1",4),"royalblue",rep("peru",2))
  ###############
  # Plot network
  XY<-list(coordX,coefY)
  myfunc<-function(x)sqrt(x)*5
  net<-plot.network(Z_cs=networks[[z]][[2]],T_cs=networks[[z]][[1]],E_cs=networks[[z]][[4]],R_cs=networks[[z]][[3]],
                    nl=2,include.detritus=TRUE, selfloop=FALSE, TP.method=1,
                    tSim=0.6, select.XY=XY, select.cols=nodcols, select.Attr=attr,labels="Default",
                    nodsize=0.25, biomass=BIOM,coefBiom=coefBIOM,coefY=NA, penwidth.func=myfunc)
  FW.plots[[z]]<-net
  #Store figure
  #   open graphical device:
  #pdf(file=paste(folders[z],seasons[z],".pdf",sep=""),         # File name
  #   width = 5, height = 6, # Width and height in inches
  #  bg = "white",          # Background color
  #  colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
  
  #render_graph(net)
  
  #   close the graphical device:
  #dev.off()
}
names(FW.plots)<-seasons

####################
#Spring 2009:
FW.plots[["SPRING2009"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[1]])
#
#Summer 2009:
FW.plots[["SUMMER2009"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[2]])
#
#Autumn 2009:
FW.plots[["AUTUMN2009"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[3]])
#
#Spring 2010:
FW.plots[["SPRING2010"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[4]])
#
#Summer 2010:
FW.plots[["SUMMER2010"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[5]])
#
#Autumn 2010:
FW.plots[["AUTUMN2010"]][["global_attrs"]][13,"value"]<-22
render_graph(FW.plots[[6]])
#