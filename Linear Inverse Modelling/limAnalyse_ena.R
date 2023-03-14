#### ENA ANALYSIS ####

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
liv<-16
nl<-2
externals<-c("Imp","Resp","Exp")

#####################################################################################
itr<-c(rep(5995,5),6994)
seasons<-c("2009 spring","2009 summer","2009 autumn","2010 spring","2010 summer","2010 autumn")
seasonfolder<-c("Spring 2009","Summer 2009","Autumn 2009","Spring 2010","Summer 2010","Autumn 2010")
folder<-"/Version 5"
all_seasons<-c("SPRING2009","SUMMER2009","AUTUMN2009","SPRING2010","SUMMER2010","AUTUMN2010")
for(z in 1:length(seasons)){
myseason<-seasons[z]
limfile<-paste(seasonfolder[z],folder,"/ENAinput/itr",itr[z],"_",all_seasons[z],"_V5.csv",sep="")
#Import LIM solutions (1000-5000 iterations):
LIM1000<-fread(file = limfile, na.strings = "", dec = "," , data.table = FALSE)

#Import information on living compartment biomass:
CBIOM<-fread(file = "Parameters/CBIOM.csv", na.strings = "", dec = "," , data.table = FALSE)
BIOM<-CBIOM[which(CBIOM[,"season"]==myseason),3:ncol(CBIOM)]
BIOM$Bac<-BIOM$Doc<-NA
#reorder biomass vector according to compartment order
BIOM<-BIOM[,comps]
ABSENTcomps<-comps[which(BIOM[1,]==0)]#check absence of compartments 
BIOM[,ABSENTcomps]<-0.001 #attribute absent compartments a minimum biomass value


#####################################################################################
#Format LIM solutions into flow matrices:
#-->stored in nested list
#each list element contains a)extended flow matrix, 
#                           b)list with diet matrix & import/export/respiration vectors
#of a certain lim solution
FROM<-c()
TO<-c()
for(i in 1:ncol(LIM1000)){
  from<-strsplit(colnames(LIM1000)[i],split="->",fixed=TRUE)[[1]][1]
  FROM<-c(FROM,from)
  to<-strsplit(colnames(LIM1000)[i],split="->",fixed=TRUE)[[1]][2]
  TO<-c(TO,to) }

LIM_ALL<-list()#list to store flow matrices
is.balanced<-rep(NA,nrow(LIM1000))#check if lim solutions are mass balanced
for(i in 1:nrow(LIM1000)){
  allflows<-unname(as.numeric(LIM1000[i,]))
  Tstar<-matrix(0,nrow=(length(comps)+3),ncol=(length(comps)+3),
               dimnames=list(c(comps,externals),c(comps,externals)))
for(u in 1:length(allflows)){
  myrow<-which(rownames(Tstar)==FROM[u])
  mycol<-which(colnames(Tstar)==TO[u])
  Tstar[myrow,mycol]<-allflows[u]}
  Tstar<-Tstar[-c((n+2),(n+3)),-c(n+1)]
  V2<-list(Tstar[1:n,1:n],Tstar[n+1,1:n],Tstar[1:n,n+1],Tstar[1:n,n+2])
  names(V2)<-c("T matrix","Z","R","E")
  bal<-sum(V2[["Z"]])-sum(V2[["R"]])-sum(V2[["E"]])
  if(round(bal,digits=6)==0){is.balanced[i]<-TRUE
  #if(bal==0){is.balanced[i]<-TRUE
  }else{is.balanced[i]<-FALSE}
  LIMsol<-list(Tstar,V2)
  names(LIMsol)<-c("Extended Transfer Matrix","Diet matrix & Inp/Resp/Exp vectors")
  LIM_ALL[[i]]<-LIMsol}

######################################################################################
#Check mass balance:
length(which(is.balanced==FALSE))
#-->inflows-outflows are balanced when rounded to the 6th digit
balanced<-"unbalanced"
########################################################################################
######### ENA CALCULATIONS

###A) DESCRIPTORS & INDICES based on unweighted flow matrix - one per season
#Connectance index: Nb of links/(Nb of species)^2
#Nb of cycles/cycle distribution

###B) DESCRIPTORS & INDICES based on weighted flow matrix: value with confidence interval

indices<-c("INTER","IMP","EXP","TSTp","TSTf","INTER/TST","IMP/TST","EXP/TST","NPP","NPP/RESP",
           "Herb","Det","Bact","APL","D/H","TE12","TE23","MTL","SOI",
           "H","Ov/DC","A/DC","R/DC","AMI","Ai/DCi","ELD","Ceff","LDq",
           "FCI","NCycles","C_cycled","C_cycled.small","C_cycled.big")

#Get functions for calculation of indices:
source("ENA PACKAGE/network_indices.R")
source("ENA PACKAGE/adapted_ENAfuncs.R")
source("ENA PACKAGE/cycling_analysis.R")

#Data frame to store indices results of all lim solutions:
INDS_ALL<-as.data.frame(matrix(NA,nrow=nrow(LIM1000),ncol=length(indices),
                               dimnames=list(NULL,indices)))#Attributes & indices
TROPH<-list()#list to store outputs from trophic analysis

#Calculate descriptors and indices:
trophic.analysis<-TRUE

for(i in 1:length(LIM_ALL)){
  Tstar<-LIM_ALL[[i]][[1]]#extended transfer matrix
  T_f1<-LIM_ALL[[i]][[2]][[1]]#matrix of internal flows
  Z_f1<-LIM_ALL[[i]][[2]][[2]]#import vector
  R_f1<-LIM_ALL[[i]][[2]][[3]]#respiration vector
  E_f1<-LIM_ALL[[i]][[2]][[4]]#export vector
  
### ATTRIBUTES ###
  
  #(1) --- Sum of internal/import/export flows
  inter<-sum(T_f1)
  imp<-sum(Z_f1)
  exp<-sum(E_f1)
  #
  #
  #(2) --- Total System Throughput and Total System Throughflow
  TSTp<-Tst(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)
  opt1<-TSTp-exp-sum(R_f1)
  opt2<-inter+imp
  TSTf<-opt1-opt2
  if(round(TSTf,digits=6)==0){TSTf<-opt1
  }else{TSTf<-NA
  print(paste(i,": Error in TSTf, opt1 unequals opt2",sep=" "))}
  #   --- Internal Throughput/TST
  inter.relative<-inter/TSTp
  #   --- Import/TST
  imp.relative<-imp/TSTp
  #   --- Export/TST
  exp.relative<-exp/TSTp
  #
  #
  #(3) --- Net Particulate & Dissolved Primary Production
  NPP<-sum(Z_f1[c("Dia","Phae")])-sum(R_f1[c("Dia","Phae")])
  #    --- NPP/Community respiration 
  resp<-NPP/sum(R_f1[c(1:liv)])
  #    --- Herbivory
  herb<-sum(T_f1[c("Dia","Phae"),1:liv])
  #    --- Detrivory
  det<-sum(T_f1[c("Doc","Poc"),1:liv])
  #    --- Bactivory
  bac<-unname(sum(T_f1["Bac",1:liv]))
  #
  #
  #(4) --- Average Path Length
  i1<-TSTf/imp
  #
  #
  #(5) --- Detrivory-Herbivory Ratio: Poc--> + Doc--> vs. Phyto--> 
  i2<-(sum(T_f1["Poc",c(1:liv)])+sum(T_f1["Doc",c(1:liv)]))/sum(T_f1[c("Dia","Phae"),c(1:liv)])
  #
  #
  #(6) --- Transfer Efficiency - require TROPHIC ANALYSIS
  if(trophic.analysis==TRUE){
  troph<-TP.CTA(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl)
  i3<-troph[["Trophic Efficiency"]][1]
  i4<-troph[["Trophic Efficiency"]][2]
  }else{i4<-i3<-NA}
  #
  #
  #(7)--- Mean Trophic Level of zooplankton/herring community
  if(trophic.analysis==TRUE){
  i5<- sum(as.numeric(troph[["Effective Trophic Levels of each species"]][3:15]*BIOM[3:15]))/sum(as.numeric(BIOM[3:15]))
  }else{i5<-NA}  
  #
  #
  #(8) --- System Omnivory Index- require TROPHIC ANALYSIS
  #Note: if TI (see function) includes zeros, SOI becomes NA because log(0) can not be computed
  #-->would not be a problem if sum(na.rm=T) is applied
  if(trophic.analysis==TRUE){
    i6<-SOI(inp=Z_f1,inter=T_f1,outt=E_f1,diss=R_f1,nl=nl)
  }else{i6<-NA}

### FLOW ORGANISATION ###
  
  #(1) --- Shannon´s Flow Diversity
  i7<-shannon.flow(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1,nl=nl,type="whole",partial.matrix=FALSE, k=1)
  #
  #
  #(2) --- Relative Overhead/Ascendency/Redundancy
  output<-Overhead(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE, show.components=TRUE)
  i8<-output[1,2]
  i9<-Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)[2]
  i10<-output[2,2]
  #
  #
  #(3) --- Average Mutual information
  i11<-Ami(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1, partial.matrix = FALSE)
  #
  #
  #(4) --- Internal Ascendency (normalized by internal development capacity)
  i12<-Internal.Asc(Z_cs=Z_f1,T_cs=T_f1,E_cs=E_f1,R_cs=R_f1)[,2]
  #
  #
  #(5) --- Effective Link Density after Ulanowicz et al. (2014)
  i13<-ELD(Z_cs=Z_f1, E_cs=E_f1, R_cs=R_f1, T_cs=T_f1)
  #
  #
  #(6) --- Effective Connectivity after Zorach & Ulanowicz (2003)
  i14<-eff.C(Z_cs=Z_f1, E_cs=E_f1, R_cs=R_f1, T_cs=T_f1)
  #
  #
  #(17) --- Link Density after Bersier et al. (2002)
  i15<-LDq(Z_cs=Z_f1, E_cs=E_f1, R_cs=R_f1, T_cs=T_f1)

### RECYCLING ###

  #(1) --- Finn Cycling Index (as fraction of TSTf)
  i16<-FCI(inp=Z_f1, inter=T_f1, outt=E_f1, diss=R_f1)[1]
  #
  #
  #(2) --- Cycling analysis
  cyc<-cycling.analysis(Z_cs=Z_f1, E_cs=E_f1, R_cs=R_f1, T_cs=T_f1, select="Default", nmax=10000)
  i17<-cyc[[1]] # number of cycles
  i18<-sum(cyc[[5]],na.rm=TRUE) # amount of matter cycled
  i19<-sum(cyc[[5]][1:2],na.rm=TRUE) # amount of matter cycled through short cycles (length = 2)
  i20<-sum(cyc[[5]][3:length(cyc[[5]])],na.rm=TRUE) # amount of matter cycled through long cycles (length > 2)

  
  
#Store results in data frame & lists:
  inds_all<-c(inter,imp,exp,TSTp,TSTf,inter.relative,imp.relative,exp.relative,NPP,resp,herb,det,bac,
              i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20)
  INDS_ALL[i,]<-inds_all
  TROPH[[i]]<-troph #trophic analysis
}

########################## CHECK NA VALUES ##############################################################
#-->TO DO!

########################## GET MEAN VALUES & 95%-CONFIDENCE INTERVALS OF ################################
#a)flows
#b)indices & descriptors
#c)trophic analysis outputs
#d)mixed trophic impact outputs
#e)interaction strength matrices
saveResults<-paste("C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Results/",seasonfolder[z],"/",sep="")
names.quantiles<-c("Q0","Q2.5","Q97.5","Q100")
foodweb.interactions<-c()
for(i in 1:length(comps)){
  IS.mycomp<-paste(comps[i],"->",comps[-c(i,length(comps))],sep=" ")
  foodweb.interactions<-c(foodweb.interactions,IS.mycomp)}

######### FLOWS #################################
FLOWS.CI95<-data.frame(Flow=colnames(LIM1000),
                       Mean=apply(LIM1000,2,mean),
                       Q0=apply(LIM1000,2,function(x)quantile(x,probs=0)),
                       Q2.5=apply(LIM1000,2,function(x)quantile(x,probs=0.025)),
                       Q97.5=apply(LIM1000,2,function(x)quantile(x,probs=0.975)),
                       Q100=apply(LIM1000,2,function(x)quantile(x,probs=1)))
#store results:
write_xlsx(FLOWS.CI95,path=paste(saveResults,"FLOWS_CI95.xlsx",sep=""))

######### INDICES & DESCRIPTORS #################
INDS.CI95<-data.frame(Index=indices,
                      Mean=apply(INDS_ALL,2,mean,na.rm=TRUE),
                      Q0=apply(INDS_ALL,2,function(x)quantile(x,probs=0,na.rm=TRUE)),
                      Q2.5=apply(INDS_ALL,2,function(x)quantile(x,probs=0.025,na.rm=TRUE)),
                      Q97.5=apply(INDS_ALL,2,function(x)quantile(x,probs=0.975,na.rm=TRUE)),
                      Q100=apply(INDS_ALL,2,function(x)quantile(x,probs=1,na.rm=TRUE)))
#store results:
write_xlsx(INDS.CI95,path=paste(saveResults,"INDS_CI95.xlsx",sep=""))

######### TROPHIC LEVELS ########################
for(i in 1:length(TROPH)){
  if(i==1){TPs<-as.data.frame(matrix(TROPH[[i]][["Effective Trophic Levels of each species"]],nrow=1,
                                     ncol=length(comps),byrow=TRUE,dimnames(list(NULL,NULL))))
  }else{
    TPs<-rbind(TPs,TROPH[[i]][["Effective Trophic Levels of each species"]])}}
TrophLEVEL.CI95<-data.frame(Comp=comps,
                            Mean=apply(TPs,2,mean,na.rm=TRUE),
                            Q0=apply(TPs,2,function(x)quantile(x,probs=0,na.rm=TRUE)),
                            Q2.5=apply(TPs,2,function(x)quantile(x,probs=0.025,na.rm=TRUE)),
                            Q97.5=apply(TPs,2,function(x)quantile(x,probs=0.975,na.rm=TRUE)),
                            Q100=apply(TPs,2,function(x)quantile(x,probs=1,na.rm=TRUE)))
#store results:
write_xlsx(TrophLEVEL.CI95,path=paste(saveResults,"TrophLEVEL_CI95.xlsx",sep=""))

######## LINDEMAN SPINE ##########################
LINDS.CI95<-vector(mode="list",length=10)
for(i in 1:length(LINDS.CI95)){
  if(i<=3){
    LINDS.CI95[[i]]<-cbind(data.frame(Mean=NA),
                           as.data.frame(matrix(NA,nrow=1, ncol=length(names.quantiles),
                                                dimnames=list(NULL,names.quantiles))))
  }else{
    LINDS.CI95[[i]]<-cbind(data.frame(TrophicLevel=paste("TL",1:20,sep=""),Mean=NA),
                           as.data.frame(matrix(NA,nrow=20, ncol=length(names.quantiles),
                                                dimnames=list(NULL,names.quantiles))))}}
names(LINDS.CI95)<-c("Detrivory", "Input to Detrital Pool", "Circulation within Detrital Pool",
                     "Canonical Exports", "Canonical Respirations", "Grazing Chain", 
                     "Returns to Detrital Pool", "Lindeman Spine", "Trophic Efficiency",
                     "Canonical Imports")#Data for Lindeman Spine
for(i in 1:length(names(LINDS.CI95))){
if(i<=3){
  for(u in 1:length(TROPH)){
    if(u==1){vec<-as.numeric(TROPH[[u]][[names(LINDS.CI95)[i]]])
    }else{vec<-c(vec,as.numeric(TROPH[[u]][[names(LINDS.CI95)[i]]]))}} 
  LINDS.CI95[[names(LINDS.CI95)[i]]][1,"Mean"]<-mean(vec,na.rm=TRUE)
  LINDS.CI95[[names(LINDS.CI95)[i]]][1,"Q0"]<-quantile(vec,probs=0,na.rm=TRUE)
  LINDS.CI95[[names(LINDS.CI95)[i]]][1,"Q2.5"]<-quantile(vec,probs=0.025,na.rm=TRUE)
  LINDS.CI95[[names(LINDS.CI95)[i]]][1,"Q97.5"]<-quantile(vec,probs=0.975,na.rm=TRUE)
  LINDS.CI95[[names(LINDS.CI95)[i]]][1,"Q100"]<-quantile(vec,probs=1,na.rm=TRUE)
}else{
  df<-as.data.frame(matrix(NA,nrow=length(TROPH),
                           ncol=20,dimnames=list(NULL,NULL),byrow=TRUE))
 for(u in 1:length(TROPH)){
   if(length(TROPH[[u]][[names(LINDS.CI95)[i]]])!=0){
   df[u,1:length(TROPH[[u]][[names(LINDS.CI95)[i]]])]<-as.numeric(TROPH[[u]][[names(LINDS.CI95)[i]]])}}
  
   LINDS.CI95[[names(LINDS.CI95)[i]]][,"Mean"]<-apply(df,2,mean,na.rm=TRUE)
   LINDS.CI95[[names(LINDS.CI95)[i]]][,"Q0"]<-apply(df,2,function(x)quantile(x,probs=0,na.rm=TRUE))
   LINDS.CI95[[names(LINDS.CI95)[i]]][,"Q2.5"]<-apply(df,2,function(x)quantile(x,probs=0.025,na.rm=TRUE))
   LINDS.CI95[[names(LINDS.CI95)[i]]][,"Q97.5"]<-apply(df,2,function(x)quantile(x,probs=0.975,na.rm=TRUE))
   LINDS.CI95[[names(LINDS.CI95)[i]]][,"Q100"]<-apply(df,2,function(x)quantile(x,probs=1,na.rm=TRUE))
}}

#store results:
write_xlsx(LINDS.CI95,path=paste(saveResults,"LINDS_CI95.xlsx",sep=""))
}#END OF SEASON LOOP
