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
  ########################################################################################
  ######### ENA CALCULATIONS: additional attributes
  
  indices<-c("GPP","BP","Egestion","Respiration","GrossHetP")
  
  #Get functions for calculation of indices:
  source("ENA PACKAGE/network_indices.R")
  source("ENA PACKAGE/adapted_ENAfuncs.R")
  source("ENA PACKAGE/cycling_analysis.R")
  
  #Data frame to store indices results of all lim solutions:
  INDS_ALL<-as.data.frame(matrix(NA,nrow=nrow(LIM1000),ncol=length(indices),
                                 dimnames=list(NULL,indices)))#Attributes & indices
  
  
  for(i in 1:length(LIM_ALL)){
    Tstar<-LIM_ALL[[i]][[1]]#extended transfer matrix
    T_f1<-LIM_ALL[[i]][[2]][[1]]#matrix of internal flows
    Z_f1<-LIM_ALL[[i]][[2]][[2]]#import vector
    R_f1<-LIM_ALL[[i]][[2]][[3]]#respiration vector
    E_f1<-LIM_ALL[[i]][[2]][[4]]#export vector
    
    ### ATTRIBUTES ###
    
    #(1) --- GPP
    gpp<-sum(Z_f1[c("Dia","Phae")])
    #
    #
    #(2) --- Bacterial Production
    bp<-sum(T_f1["Bac",])-T_f1["Bac","Poc"]+E_f1["Bac"]
    #
    #
    #(3) --- Sum of egestions
    poc<-sum(T_f1[,"Poc"])
    #
    #
    #(4) --- Respiration
    resp<-sum(R_f1)
    #
    #
    #(5) --- Gross Heterotroph Production
    GhetP<-sum(T_f1[,living[3:length(living)]])+sum(Z_f1[living[3:length(living)]])
    
    
    #Store results in data frame & lists:
    inds_all<-c(gpp,bp,poc,resp,GhetP)
    INDS_ALL[i,]<-inds_all
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
  
  ######### additional ATTRIBUTES #################
  INDS.CI95<-data.frame(Index=indices,
                        Mean=apply(INDS_ALL,2,mean,na.rm=TRUE),
                        Q0=apply(INDS_ALL,2,function(x)quantile(x,probs=0,na.rm=TRUE)),
                        Q2.5=apply(INDS_ALL,2,function(x)quantile(x,probs=0.025,na.rm=TRUE)),
                        Q97.5=apply(INDS_ALL,2,function(x)quantile(x,probs=0.975,na.rm=TRUE)),
                        Q100=apply(INDS_ALL,2,function(x)quantile(x,probs=1,na.rm=TRUE)))
  #store results:
  write_xlsx(INDS.CI95,path=paste(saveResults,"ATTR_CI95.xlsx",sep=""))
}#END OF SEASON LOOP  
