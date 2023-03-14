#### EXPLORE LIM MODELS AND NETWORK FLOWS ####

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/")

#load packages:
library("data.table")
library("ggplot2")
library("writexl")
library("LIM")

i<-6
seasonfolder<-c("Spring 2009","Summer 2009","Autumn 2009","Spring 2010","Summer 2010","Autumn 2010")
all_seasons<-c("SPRING2009","SUMMER2009","AUTUMN2009","SPRING2010","SUMMER2010","AUTUMN2010")
myseason<-seasonfolder[i]
folder<-"/Version 5"                 
LIMinput<-paste(all_seasons[i],"_version5.input",sep="")
LIMfile1<-paste("itr1000",all_seasons[i],"_V5central.csv",sep="")
LIMfile2<-paste("itr1000",all_seasons[i],"_V5lsei.csv",sep="")
LIMfile3<-paste("itr1000",all_seasons[i],"_V5specific1.csv",sep="")
LIMfile4<-paste("itr1000",all_seasons[i],"_V5specific2.csv",sep="")
LIMfile5<-paste("itr1000",all_seasons[i],"_V5specific3.csv",sep="")
LIMfile6<-paste("itr1000",all_seasons[i],"_V5specific4.csv",sep="")
LIMfile7<-paste("itr1000",all_seasons[i],"_V5specific5.csv",sep="")
LIMfile8<-NULL#paste("itr1000",all_seasons[i],"_V5specific6.csv",sep="")

#Import LIM solutions (1000-5000 iterations):
LIM1000central<-fread(file = paste(myseason,folder,"/",LIMfile1,sep=""), 
               na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-LIM1000central
if(length(LIMfile2)!=0){
LIM1000lsei<-fread(file = paste(myseason,folder,"/",LIMfile2,sep=""), 
                   na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),])}
if(length(LIMfile3)!=0){
LIM1000spec1<-fread(file = paste(myseason,folder,"/",LIMfile3,sep=""), 
                   na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),])}
if(length(LIMfile4)!=0){
LIM1000spec2<-fread(file = paste(myseason,folder,"/",LIMfile4,sep=""), 
                   na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),],LIM1000spec2[2:nrow(LIM1000central),])}
if(length(LIMfile5)!=0){
LIM1000spec3<-fread(file = paste(myseason,folder,"/",LIMfile5,sep=""), 
                   na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),],LIM1000spec2[2:nrow(LIM1000central),],LIM1000spec3[2:nrow(LIM1000central),])}
if(length(LIMfile6)!=0){
LIM1000spec4<-fread(file = paste(myseason,folder,"/",LIMfile6,sep=""), 
                      na.strings = "", dec = "," , data.table = FALSE)
LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),],LIM1000spec2[2:nrow(LIM1000central),],LIM1000spec3[2:nrow(LIM1000central),],LIM1000spec4[2:nrow(LIM1000central),])}
if(length(LIMfile7)!=0){
  LIM1000spec5<-fread(file = paste(myseason,folder,"/",LIMfile7,sep=""), 
                      na.strings = "", dec = "," , data.table = FALSE)
  LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),],LIM1000spec2[2:nrow(LIM1000central),],
                 LIM1000spec3[2:nrow(LIM1000central),],LIM1000spec4[2:nrow(LIM1000central),],LIM1000spec5[2:nrow(LIM1000central),])}
if(length(LIMfile8)!=0){
  LIM1000spec6<-fread(file = paste(myseason,folder,"/",LIMfile8,sep=""), 
                      na.strings = "", dec = "," , data.table = FALSE)
  LIM1000<-rbind(LIM1000central,LIM1000lsei[2:nrow(LIM1000central),],LIM1000spec1[2:nrow(LIM1000central),],LIM1000spec2[2:nrow(LIM1000central),],
                 LIM1000spec3[2:nrow(LIM1000central),],LIM1000spec4[2:nrow(LIM1000central),],LIM1000spec5[2:nrow(LIM1000central),],
                 LIM1000spec6[2:nrow(LIM1000central),])}

#Save model output:
ENAinput<-paste(all_seasons[i],"_V5.xlsx",sep="")
write_xlsx(LIM1000,path=paste(myseason,folder,"/ENAinput/itr",nrow(LIM1000),"_",ENAinput,sep=""))
###################################################################################
#####CHECK 1: COVERAGE OF SOLUTION SPACE

#check coverage of the solution space (cf. jumpsize):
readLIM<-Read(paste(myseason,folder,"/",LIMinput,sep=""))
LIM<-Setup(readLIM) 
FlowRanges<-xranges(E=LIM$A,F=LIM$B,G=LIM$G,H=LIM$H,
                    ispos=TRUE, central = TRUE)
rownames(FlowRanges)<-LIM$Unknowns
SampleRange <- cbind(data.frame(FlowNb=1:length(LIM$Unknowns)),FlowID=LIM$Unknowns,data.frame(FlowRanges))
SampleRange$samplemin <- apply(LIM1000, 2, min)
SampleRange$samplemax <- apply(LIM1000, 2, max)
SampleRange$percCover <- ((SampleRange$samplemax - SampleRange$samplemin) / 
                            (SampleRange$max - SampleRange$min) * 100)
SampleLOWS<-SampleRange[which(SampleRange$percCover<50),]

# Mean percentage covered range
mean(SampleRange$percCover, na.rm = T)
# % of flows with covered range >75%
length(which(SampleRange$percCover>=75))/length(LIM$Unknowns)

###############################################################################################
### CALCULATE MEAN + SD VALUE OF EACH FLOW + 95% CONFIDENCE INTERVALS
meanvalues <- cbind(LIM$Unknowns, colMeans(LIM1000))
standarddev <- cbind(LIM$Unknowns, sqrt(diag(var(LIM1000))))
Summary <- data.frame(Flow=LIM$Unknowns, MEAN=as.numeric(meanvalues[,2]), 
                      SD=as.numeric(standarddev[,2]))
rownames(Summary)<-1:length(LIM$Unknowns)
Summary$LowerCI95<-apply(LIM1000,2,FUN=function(x)quantile(x,probs=0.025))
Summary$UpperCI95<-apply(LIM1000,2,FUN=function(x)quantile(x,probs=0.975))
Summary<-cbind(Summary, SampleRange[,c("min","max","samplemin","samplemax","percCover")])
saveOutput<-paste("C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/",myseason,folder,"/Summary_itr",nrow(LIM1000),".xlsx",sep="")
write_xlsx(Summary[order(Summary[,"percCover"],decreasing=TRUE),], path=saveOutput)
Summary<-Summary[order(Summary[,"percCover"],decreasing=TRUE),]
####################################################################################
#### CHECK 2: CONVERGENCE
#-->visually checked for 10 randomly chosen flows

# Step size used for sampling all random samples. Guideline:
# If iterations >1000 use iterations / 1000.
# Must be integer.
step <- nrow(LIM1000)/1000
itr <- nrow(LIM1000)
checkFlows<-sample(x=1:147, size=10, replace=FALSE)
for(i in 1:length(checkFlows)){
# The nr of the flow you want to check.
flownr  <- checkFlows[i]

###################
# Check Convergence
#for a given flow
###################

nsample <- itr/step
sdevs   <- vector("numeric", nsample)
means   <- vector("numeric", nsample)
for (i in 1:nsample){
  random   <- sort(runif(itr), index.return=TRUE)
  sdevs[i] <-   sd(LIM1000[random$ix[1:(i*step)], flownr])
  means[i] <- mean(LIM1000[random$ix[1:(i*step)], flownr])
}
plot(sdevs)
plot(means)

###################################
# Calculate fluctuation reduction
# for a given flow
###################################

# Error margin of 2% of average stand. dev.
stdevofflow <- as.numeric(standarddev[flownr, 2])
errormargin <- stdevofflow * 0.02
plot(sdevs)
abline(h = stdevofflow + errormargin/2)
abline(h = stdevofflow - errormargin/2)
}

####################################################################################
#### CHECK 3: ENERGY BUDGETS OF COMPARTMENTS

#Compartments:
comps<-c("Bac","Dia","Phae","Cil","Nsci","Tun","Clado","Cop","Biv","Gastr","Poly","Hydro","Ppil","Mlei","Bcu","Her","Doc","Poc")
funcGroups<-c(rep("bottom",3),"mic",rep("meso",4),rep("mero",3),rep("top",5),rep("nl",2))

#Create vectors to identify incoming/outgoing compartments:
FROM<-c()
TO<-c()
for(i in 1:ncol(LIM1000)){
  from<-strsplit(colnames(LIM1000)[i],split="->",fixed=TRUE)[[1]][1]
  FROM<-c(FROM,from)
  to<-strsplit(colnames(LIM1000)[i],split="->",fixed=TRUE)[[1]][2]
  TO<-c(TO,to) }
#Create data frame indicating AVERAGE + SD import/specific ingestion/total ingestion/
#total production/export/egestion/respiration:
EBUDGETS.mean<-EBUDGETS.sd<-
  as.data.frame(matrix(NA,nrow=length(comps),ncol=(length(comps)+9),
      dimnames=list(NULL,c("Comp","funcGroup","Imp",comps,"Ingestion","Production",
                           "Export","DOCexs","Egestion","Respiration"))))
EBUDGETS.mean$Comp<-EBUDGETS.sd$Comp<-comps
EBUDGETS.mean$funcGroup<-EBUDGETS.sd$funcGroup<-funcGroups
for(i in 1:length(comps)){
  IN<-which(TO==comps[i])
  OUT<-which(FROM==comps[i])
  LOSS<-c(which(TO[OUT]=="Doc"),which(TO[OUT]=="Poc"),which(TO[OUT]=="Resp"))
  ebudget.mean<-apply(LIM1000[,c(IN,OUT)],2,mean)
  ebudget.sd<-apply(LIM1000[,c(IN,OUT)],2,sd)
  ebudget<-as.data.frame(rbind(FROM[c(IN,OUT)],TO[c(IN,OUT)],ebudget.mean,ebudget.sd))
  rownames(ebudget)<-c("FROM","TO","MEAN","SD")
  ING<-LIM1000[,IN]
  if(class(ING)=="data.frame")ING<-apply(ING,1,sum)
  ingestion.mean<-mean(ING)
  ingestion.sd<-sd(ING)
  PROD<-LIM1000[,OUT]
  if(length(LOSS)>0)PROD<-PROD[,-LOSS]
  if(class(PROD)=="data.frame")PROD<-apply(PROD,1,sum)
  production.mean<-mean(PROD)
  production.sd<-sd(PROD)
diet<-as.character(unname(ebudget["FROM",which(ebudget["TO",]==comps[i])]))  
EBUDGETS.mean[i,diet]<-ebudget["MEAN",which(ebudget["TO",]==comps[i])]
EBUDGETS.sd[i,diet]<-ebudget["SD",which(ebudget["TO",]==comps[i])]  
if(funcGroups[i]!="nl"){
EBUDGETS.mean[i,"Ingestion"]<-ingestion.mean
EBUDGETS.sd[i,"Ingestion"]<-ingestion.sd}
EBUDGETS.mean[i,"Production"]<-production.mean
EBUDGETS.sd[i,"Production"]<-production.sd
if(length(which(ebudget["TO",]=="Exp"))>0){
EBUDGETS.mean[i,"Export"]<-ebudget["MEAN",which(ebudget["TO",]=="Exp")]
EBUDGETS.sd[i,"Export"]<-ebudget["SD",which(ebudget["TO",]=="Exp")]}
if(length(which(ebudget["TO",]=="Doc"))>0&& comps[i]!="Doc"){
  EBUDGETS.mean[i,"DOCexs"]<-ebudget["MEAN",which(ebudget["TO",]=="Doc")]
  EBUDGETS.sd[i,"DOCexs"]<-ebudget["SD",which(ebudget["TO",]=="Doc")]}
if(length(which(ebudget["TO",]=="Poc"))>0 && comps[i]!="Poc"){
EBUDGETS.mean[i,"Egestion"]<-ebudget["MEAN",which(ebudget["TO",]=="Poc")]
EBUDGETS.sd[i,"Egestion"]<-ebudget["SD",which(ebudget["TO",]=="Poc")]}
if(length(which(ebudget["TO",]=="Resp"))>0){
EBUDGETS.mean[i,"Respiration"]<-ebudget["MEAN",which(ebudget["TO",]=="Resp")]
EBUDGETS.sd[i,"Respiration"]<-ebudget["SD",which(ebudget["TO",]=="Resp")]}
EBUDGETS.mean[i,which(is.na(EBUDGETS.mean[i,])==TRUE)]<-0
EBUDGETS.sd[i,which(is.na(EBUDGETS.sd[i,])==TRUE)]<-0}
for(i in 3:ncol(EBUDGETS.mean)){
  EBUDGETS.mean[,i]<-as.numeric(EBUDGETS.mean[,i])
  EBUDGETS.sd[,i]<-as.numeric(EBUDGETS.sd[,i])}

##############################################################
### PLOT ENERGY BUDGETS ###
limPLOT<-list()
for(i in 1:length(unique(funcGroups))){
  MEANs<-as.data.frame(EBUDGETS.mean[which(EBUDGETS.mean[,"funcGroup"]==unique(funcGroups)[i]),])
  SDs<-as.data.frame(EBUDGETS.sd[which(EBUDGETS.sd[,"funcGroup"]==unique(funcGroups)[i]),])
if(funcGroups[i]!="nl"){
  outs<-unname(which(apply(MEANs[,3:ncol(MEANs)],2,sum)==0))
  FlowTypes<-colnames(MEANs)[3:ncol(MEANs)][-outs]
}else{FlowTypes<-c("Imp","Production","Export")}
plotEB<-as.data.frame(matrix(NA,nrow=0,ncol=4,dimnames=list(NULL, c("Comp",
                                          "FlowType","Mean","Sd"))))
for(u in 1:nrow(MEANs)){
  flux<-data.frame(Comp=rep(MEANs[u,"Comp"],length(FlowTypes)),FlowType=FlowTypes,
                   Mean=NA,Sd=NA)
  for(s in 1:length(FlowTypes)){
    flux[s,"Mean"]<-MEANs[u,FlowTypes[s]]
    flux[s,"Sd"]<-SDs[u,FlowTypes[s]]}
  plotEB<-rbind(plotEB,flux)}
#plotEB$Fmin<-plotEB$Mean-plotEB$Sd
#plotEB$Fmax<-plotEB$Mean+plotEB$Sd
if(funcGroups[i]!="nl"){
  cols<-rep("slateblue2",length(FlowTypes))
  if(length(which(FlowTypes=="DOCexs"))==0){
    cols[(length(FlowTypes)-3):length(FlowTypes)]<-rep("tomato2",4)
  }else{cols[(length(FlowTypes)-4):length(FlowTypes)]<-rep("tomato2",5)}
}else{cols<-rep("darkgoldenrod",length(Flowtypes))}
plotEB$FlowType<-factor(plotEB$FlowType,levels=FlowTypes)
if(length(unique(MEANs$Comp))>1){
plot<-
  ggplot(plotEB,aes(x=FlowType,y=Mean))+
  theme_bw()+
  geom_bar(aes(fill=FlowType),position="dodge", stat = "identity",show.legend=FALSE) +
  scale_fill_manual(values=cols,breaks=FlowTypes)+
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), color="black", width=.2,
                position=position_dodge(.9))+
  ylab("mg C/m^3/day")+
  theme(panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line( size=.1, color="grey" ),
        legend.title = element_blank(),
        axis.text.x = element_text(face="bold",angle=30 ,size=9, color="grey48"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10, color="grey28"),
        axis.title.y = element_text(face="bold", size=10, color="grey48"))+
    facet_wrap(~Comp,ncol=1,scales="free")
}else{
  plot<-
    ggplot(plotEB,aes(x=FlowType,y=Mean))+
    theme_bw()+
    geom_bar(aes(fill=FlowType),position="dodge", stat = "identity",show.legend=FALSE) +
    scale_fill_manual(values=cols,breaks=FlowTypes)+
    geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), color="black", width=.2,
                  position=position_dodge(.9))+
    ylab("mg C/m^3/day")+
  theme(panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line( size=.1, color="grey" ),
        legend.title = element_blank(),
        axis.text.x = element_text(face="bold",angle=30 ,size=9, color="grey48"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", size=10, color="grey28"),
        axis.title.y = element_text(face="bold", size=10, color="grey48")) }
 limPLOT[[i]]<-plot}
#Store plots as pdf files:
plotname<-c("FoodWebBase","Ciliates","Mesozooplankton","Meroplankton","TopPredator","NonLiving")
for(i in 1:length(limPLOT)){
#   open graphical device:
  pdf(file=paste(myseason,folder,"/itr",nrow(LIM1000),"_",plotname[i],".pdf",sep=""),         # File name
      width = 5, height = 6, # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
  
plot(limPLOT[[i]])
  
#   close the graphical device:
  dev.off() 
}

for(i in 1:length(limPLOT)){
  plot(limPLOT[[i]])
}
