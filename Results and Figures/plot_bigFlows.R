##### ADDITIONAL PLOTS #####

#set working directory
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/ENA analysis/Results/")

#load packages:
library("data.table")
library(doBy)
library("ggplot2")

#Seasons:
seasons<-c("spring 2009","summer 2009","autumn 2009","spring 2010","summer 2010","autumn 2010")

#Import dataframe with the seasonally 5 largest flows:
BIGs<-fread(file = "bigFlows.csv", na.strings = "", dec = "," , data.table = FALSE)

#Convert flows as % of TST:
BIGs$pct_TST<-BIGs$Mean/BIGs$TST_mean*100
#Calculate missing percentage of TST:
#Sum<-summaryBy(formula=pct_TST~Season, data=BIGs, FUN=sum)
#Sum$pct_left<-(Sum$pct_TST.sum-1)*(-1)
#Sum<-data.frame(Season=Sum$Season,Flow=rep("other flows",6),pct_TST=Sum$pct_left)

#Define order of flows:
unique(BIGs$Flow)
#Flow.order<-c("Imp->Poc", "Poc->Exp", "Imp->Phae", "Phae->Exp", "Imp->Dia",  "Dia->Exp", "Dia->Resp",  "Imp->Nsci", "Nsci->Poc", 
            #  "Doc->Bac")
Flow.order<-c("Imp->Poc", "Poc->Exp", "Imp->Phae", "Phae->Exp", "Imp->Dia",  "Dia->Exp",  "Imp->Nsci", "Nsci->Poc")
Flow.order<-Flow.order[8:1]

#Define color of flows:
#cols<-c("darkgoldenrod","khaki3","wheat4","snow3","darkcyan","skyblue","paleturquoise","yellowgreen","yellow3",
      #  "mediumorchid")
cols<-c("darkgoldenrod","khaki3","wheat4","snow3","darkcyan","skyblue","yellowgreen","yellow3")
cols<-cols[8:1]
#Plot big flows as % of TST:
setwd("C:/Hannah/Biological Oceanography/Master Thesis/Project/RESULTS/Figures/")
BIGs$Season<-factor(BIGs$Season,levels=seasons)
BIGs$Flow<-factor(BIGs$Flow,levels=Flow.order)

pdf(file="BigFlows.pdf",         # File name
    width = 6, height = 4.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(BIGs, aes(fill=Flow, y=pct_TST, x=Season))+ 
  theme_minimal()+
  geom_bar(position = "stack", stat="identity", width = 0.45, color="black")+
  scale_fill_manual(values=cols,
                    breaks=Flow.order,
                    labels=Flow.order)+
  ylim(0,100)+
  ylab("% TST")+
  xlab("")+
  theme(legend.title=element_text(size=11, face="bold", color="gray21"),
        axis.title.y=element_text(size=10, face="bold", color="black"),
        axis.text.x=element_text(angle=25,size=9, color="black"),
        panel.grid.major = element_line(colour="lightgrey", linetype="dashed"), 
        panel.grid.minor = element_blank())
# close the graphical device:
dev.off()

#POC export and import across seasons
POC<-BIGs[which(BIGs$Flow=="Imp->Poc"),]
POC<-rbind(POC, BIGs[which(BIGs$Flow=="Poc->Exp"),])
POC$Flow<-factor(POC$Flow, levels=c("Imp->Poc","Poc->Exp"))
POC$Season<-factor(POC$Season, levels=seasons)
pdf(file="POC.pdf",         # File name
    width = 6, height = 4.5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" )   # Color model (cmyk is required for most publications)
ggplot(POC, aes(fill=Flow, y=Mean, x=Season,ymin=Q2.5, ymax=Q97.5))+ 
  theme_minimal()+
  geom_bar(aes(fill=Flow),position = "dodge", stat="identity", width = 0.45, color="black")+
  geom_errorbar(color="black",position=position_dodge(0.5), width=0.2)+
  scale_fill_manual(values=cols[7:8],
                    breaks=Flow.order[7:8],
                    labels=Flow.order[7:8])+
  ylab(expression(paste("mg C/m" ^ "3","/day",sep="")))+
  xlab("")+
  theme(legend.title=element_text(size=11, face="bold", color="gray21"),
        axis.title.y=element_text(size=10, face="bold", color="black"),
        axis.text.x=element_text(angle=25,size=9, color="black"),
        panel.grid.major = element_line(colour="lightgrey", linetype="dashed"), 
        panel.grid.minor = element_blank())
# close the graphical device:
dev.off()
