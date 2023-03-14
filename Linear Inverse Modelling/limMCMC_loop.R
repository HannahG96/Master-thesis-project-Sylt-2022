######################################
### LIM MODELLING ####################
### MCMC SAMPLING - LOOP ####################

#load packages:
library("LIM")
library("limSolve")
library("writexl")
library("data.table")

working.dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/"
setwd(working.dir)

itr<-1000
seasons<-c("SPRING2009","SUMMER2009","AUTUMN2009","SPRING2010","SUMMER2010","AUTUMN2010")
Xzeros<-c("central","lsei","ldei")

for(i in c(2,3,5)){
  myseason<-seasons[i]

#Import the lim file:
readLIM<-Read(paste(myseason,"_version5.input",sep=""))
LIM<-Setup(readLIM)   

# Find the solution range of each flow:
FlowRanges<-xranges(E=LIM$A,F=LIM$B,G=LIM$G,H=LIM$H,
                    ispos=TRUE, central = TRUE)
rownames(FlowRanges)<-LIM$Unknowns

#Perform MCMC sampling for 3 different initial solutions:
for(u in 1:length(Xzeros)){
  if(Xzeros[u]=="central"){
    myx0<-unname(FlowRanges[,3])
  }else if(Xzeros[u]=="lsei"){
    myx0<-lsei(E=LIM$A, F=LIM$B,G=LIM$G,H=LIM$H,
               tol = sqrt(.Machine$double.eps), verbose = TRUE)$X     
  }else if(Xzeros[u]=="ldei"){
    myx0<-ldei(E=LIM$A, F=LIM$B,G=LIM$G,H=LIM$H,
                  tol = sqrt(.Machine$double.eps), verbose = TRUE)$X}
jumpsize<-100
jump<-unname((FlowRanges[,2] - FlowRanges[,1]))/jumpsize

XS_LIM<-xsample(E=LIM$A, F=LIM$B, A=NULL,
                  B=NULL,jmp=jump,G=LIM$G,H=LIM$H, x0=myx0,
                  iter=itr,type="mirror",outputlength=itr, test=TRUE, burninlength=1)
xs_LIM<-as.data.frame(XS_LIM$X)
colnames(xs_LIM)<-LIM$Unknowns

#Store results as excelfile:
outputname<-paste("itr",itr,myseason,"_V5",Xzeros[u],".xlsx",sep="")
saveOutput<-paste(working.dir,outputname,sep="")
write_xlsx(xs_LIM, path=saveOutput)}
}