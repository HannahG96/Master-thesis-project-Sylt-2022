################################################################
#### CREATE SPECIFIC INITIAL SOLUTIONS FOR MCMC SAMPLING #######
################################################################

#load packages:
library("LIM")
library("limSolve")
library("writexl")
library("data.table")

working.dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/R Code/"

#set working directory:
setwd(working.dir)

i<-6
folder<-"/Version 5"  
seasonfolder<-c("Spring 2009","Summer 2009","Autumn 2009","Spring 2010","Summer 2010","Autumn 2010")
all_seasons<-c("SPRING2009","SUMMER2009","AUTUMN2009","SPRING2010","SUMMER2010","AUTUMN2010")
jumps<-c(500,100,100,100,100,500)

#for(i in 6){
myseason<-seasonfolder[i]
LIMinput<-paste(myseason,folder,"/",all_seasons[i],"_version5.input",sep="")

#Parameters for calculation of specific initial solution:
u<-5
limX0<-c(paste(myseason,folder,"/Xzeros/",all_seasons[i],"_version5_x01.input",sep=""),
         paste(myseason,folder,"/Xzeros/",all_seasons[i],"_version5_x02.input",sep=""),
         paste(myseason,folder,"/Xzeros/",all_seasons[i],"_version5_x03.input",sep=""),
         paste(myseason,folder,"/Xzeros/",all_seasons[i],"_version5_x04.input",sep=""),
         paste(myseason,folder,"/Xzeros/",all_seasons[i],"_version5_x05.input",sep=""))
Xzeros<-c("specific1","specific2","specific3","specific4", "specific5")
#Import the lim file:
readLIM<-Read(LIMinput)
LIM<-Setup(readLIM) 
# Find the solution range of each flow:
FlowRanges<-xranges(E=LIM$A,F=LIM$B,G=LIM$G,H=LIM$H,
                    ispos=TRUE, central = TRUE)
rownames(FlowRanges)<-LIM$Unknowns

#Define specific initial solution to reach upper extremes of low sampled flows:
#Focus: (1) Poc -> Exp
#       (2) Nsci -> Resp
#       (3) Bac -> Resp (?)
#LOWS<-fread(file = paste(myseason,folder,"/Summary_itr1999.csv",sep=""), 
   #         na.strings = "", dec = "," , data.table = FALSE)
LOWS<-fread(file = paste(myseason,folder,"/Summary_itr5995.csv",sep=""), 
         na.strings = "", dec = "," , data.table = FALSE)
LOWS<-LOWS[which(LOWS[,"percCover"]<60),c("Flow","min","max","samplemin","samplemax")]
LOWS$x0_95pct<-LOWS$min + (LOWS$max-LOWS$min) * 0.95
LOWS$x0_75pct<-LOWS$min + (LOWS$max-LOWS$min) * 0.75
LOWS$x0_45pct<-LOWS$min + (LOWS$max-LOWS$min) * 0.45 
LOWS$x0_25pct<-LOWS$min + (LOWS$max-LOWS$min) * 0.25 
LOWS$x0_5pct<-LOWS$min + (LOWS$max-LOWS$min) * 0.05 
print(LOWS)

for(u in 5){
readLIM.x0<-Read(limX0[u])
LIM.x0<-Setup(readLIM.x0) 
#using central value of flow ranges:
x0.specific1<-xranges(E=LIM.x0$A,F=LIM.x0$B,G=LIM.x0$G,H=LIM.x0$H,ispos=TRUE, central = TRUE)[,3]
xranges(E=LIM$A,F=LIM$B,G=LIM$G,H=LIM$H,ispos=TRUE, central = TRUE)[,3]-x0.specific1
#using lsei:
x0.specific2<-lsei(E=LIM.x0$A, F=LIM.x0$B,G=LIM.x0$G,H=LIM.x0$H,
                   tol = sqrt(.Machine$double.eps), verbose = TRUE)$X
#xranges(E=LIM$A,F=LIM$B,G=LIM$G,H=LIM$H,ispos=TRUE, central = TRUE)[,3]-x0.specific2
#->DECIDING FOR LSEI:
myx0<-x0.specific2

##Sample the solution space:
itr<-1000
jumpsize<-jumps[i]
jump<-unname((FlowRanges[,2] - FlowRanges[,1]))/jumpsize

XS_LIM<-xsample(E=LIM$A, F=LIM$B, A=NULL,
                B=NULL,jmp=jump,G=LIM$G,H=LIM$H, x0=myx0,
                iter=itr,type="mirror",outputlength=itr, test=TRUE, burninlength=1)
xs_LIM<-as.data.frame(XS_LIM$X)
colnames(xs_LIM)<-LIM$Unknowns

#Store results as excelfile:
outputname<-paste(myseason,folder,"/itr",itr,all_seasons[i],"_V5",Xzeros[u],".xlsx",sep="")
saveOutput<-paste(working.dir,outputname,sep="")
write_xlsx(xs_LIM, path=saveOutput)
}#limX0 loop
#}#season loop
