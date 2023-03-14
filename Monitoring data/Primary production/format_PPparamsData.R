#######################################################
#ANALYSE & FORMAT PRIMARY PRODUCTION TIME SERIES OF SRB 
#######################################################
#-->sampled at Königshafen

#set working directory
working_dir<-"C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/"
setwd(working_dir)

#load packages
library(ggplot2)
library(writexl)

#Import data (stored as tab files):

source("R functions/import_pangaea_file.R")#function to convert tab file into data frame
folder<-"PP_params_R.&H.Asmus/"
a<-rep("List_Koenigshafen_phytoplank_PP_",8)
timeseries<-c(2007:2014)
filenames<-vector()#Create vector of the filenames to be imported
for(i in 1:length(timeseries)){
  myfilename<-paste(folder,a[i],timeseries[i],".tab",sep="")
  filenames<-c(filenames,myfilename)}

datall.list<-list()#list to store all imported data files
datall.names<-vector()#vector to attribute a name to each data file
datall.colnames<-vector()#column names of all data files
stations<-rep("List Koenigshafen",8)
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
datall.colnames<-unique(datall.colnames)#unique column names of all data files
colnames<-datall.colnames[c(12,11,1:10)]#reorder column names
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

#Format PP parameter columns into numeric
for(i in 5:length(colnames))datall[,i]<-as.numeric(datall[,i])
####Note: in PPparams tab file 2012/13 --> missing values, NAs introduced

#Add month column
datall$month<-NA
for(i in 1:nrow(datall))datall[i,"month"]<-strsplit(as.character(datall[i,"Date_Time"]),split="-")[[1]][2]
datall$month<-as.numeric(datall$month)

#Unit Conversion of relevant parameters:
#-->Biomass (of tot phytoplankton or tot plankton?????????) 
#-->GPP [mg/l/h] in mg/m3/day
datall$GPP_C<-datall$GPP_C*1000*24
#-->NPP [mg POC+DOC/l/h] in mg POC+DOC/m3/day
datall$NP_DOC_POC<-datall$NP_DOC_POC*1000*24
#-->Respiration [umol O2/l/h] in mg C/m3/day
#after Brey (2010): mmol O2-->mg O2 = 33.191 ; mg O2-->mg C = 0.309
datall$Resp_O2<-(datall$Resp_O2*33.191*0.309*24)*(-1)

###################################################################################
###CALCULATE PARAMETER ESTIMATES
#What do I need? -->Net & Gross Primary Production from spring 2009-winter 2010/11
source("R functions/get_estimates_MULTIparam.R")
PPparams<-c("Biom_C","Resp_O2","GPP_C","NP_DOC_POC")
data.PPparams<-get_estimates_MULTIparam(dat=datall,timeseries=2009:2011,
                           stations=c("KH"),colnames.params=PPparams)[[1]]
#data.PPparams<-data.PPparams[1:26,] #remove March-Dec 2011 observations

#Explore seasonal patterns of primary production and respiration in SRB
df<-data.frame(year=rep(as.factor(data.PPparams$year),4),month=rep(data.PPparams$month,4),
               param=c(rep("Biomass [mg C/m3]",nrow(data.PPparams)),rep("GPP [mg C/m3/day]",nrow(data.PPparams)),
        rep("NPP of POC+DOC [mg POC+DOC/m3/day]",nrow(data.PPparams)),rep("Respiration [mg C/m3/day]",nrow(data.PPparams))), 
               rate=c(data.PPparams[,4],data.PPparams[,6],data.PPparams[,7],data.PPparams[,5]))
df$month<-factor(df$month,levels=1:12)
ggplot(df,aes(x=month,y=rate,color=year,group=year))+
  theme_classic()+
  geom_line()+
  #geom_smooth(mapping=aes(y=rate,x=year,color=season),data=df,method="loess",size=0.5, linetype="dashed",se=FALSE)+
  scale_color_manual("", values=c("mediumaquamarine","lawngreen","maroon"),
                     breaks=c("2009","2010","2011"))+
  #scale_shape_manual("", values=c(15,16,17,18,3),
                     #breaks=c("winter","spring","summer","autumn","annual"))+
  ylab("")+
  theme(plot.title=element_text(face="bold",size=10,color="navy",hjust=0.5),
        axis.text.x=element_text(angle=45,size=8,face="bold",color="slategray3"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=9,color="slategray3"))+
  facet_wrap(~param,ncol=1, scales="free")

######################################################################
#CALCULATE ESTIMATES OF SEASONAL GPP & NPP: spring 2009-winter 2010/11 
#FOR SETTING UPPER LOWER BOUNDARIES OF PHYTOPLANKTON PARAMETERS:

#monthly average GPP and GPP-NPP:
GPP_NPP2009_11<-datall[which(datall$Date_Time>="2009-03-09"),]
GPP_NPP2009_11<-GPP_NPP2009_11[which(GPP_NPP2009_11$Date_Time<="2011-02-15"),]
GPP_NPP2009_11$GPPminusNPP<-GPP_NPP2009_11$GPP_C-GPP_NPP2009_11$NP_DOC_POC
GPP_NPP2009_11<-get_estimates_MULTIparam(dat=GPP_NPP2009_11,timeseries=2009:2011,
                                        stations=c("KH"),colnames.params=c("GPP_C","NP_DOC_POC","GPPminusNPP"))[[1]]
colnames(GPP_NPP2009_11)[which(colnames(GPP_NPP2009_11)=="GPPminusNPP_monthly_AVG")]<-"GPPminusNPP"
GPP_NPP2009_11$season<-rep(c(rep("winter",2),rep("spring",3),rep("summer",3),rep("autumn",3),"winter"),3)
GPP_NPP2009_11$yearS<-GPP_NPP2009_11$year
GPP_NPP2009_11[which(GPP_NPP2009_11$month==12),"yearS"]<-GPP_NPP2009_11[which(GPP_NPP2009_11$month==12),"yearS"]+1
GPP_NPP2009_11$season<-paste(GPP_NPP2009_11$yearS,GPP_NPP2009_11$season,sep=" ")
#seasonal mins and max based on monthly averages:
GPP_NPP<-data.frame(year=c(rep(2009,3),rep(2010,4),2011),
                    season=paste(c(rep(2009,3),rep(2010,4),2011),c("spring","summer","autumn","winter"),sep=" "),
                    maxGPP=NA,minGPP=NA,meanGPP=NA,maxNPP=NA,minNPP=NA,meanNPP=NA,
                    maxGPPminusNPP=NA,minGPPminusNPP=NA,meanGPPminusNPP=NA)

for(i in 1:nrow(GPP_NPP)){
  myseason<-GPP_NPP[i,"season"]
  myvals<-GPP_NPP2009_11[which(GPP_NPP2009_11$season==myseason),]
  GPP_NPP[i,"maxGPP"]<-max(myvals$GPP_C_monthly_AVG,na.rm=TRUE)[1]
  GPP_NPP[i,"minGPP"]<-min(myvals$GPP_C_monthly_AVG,na.rm=TRUE)[1]
  GPP_NPP[i,"meanGPP"]<-mean(myvals$GPP_C_monthly_AVG,na.rm=TRUE)
  GPP_NPP[i,"maxNPP"]<-max(myvals$NP_DOC_POC_monthly_AVG,na.rm=TRUE)[1]
  GPP_NPP[i,"minNPP"]<-min(myvals$NP_DOC_POC_monthly_AVG,na.rm=TRUE)[1]
  GPP_NPP[i,"meanNPP"]<-mean(myvals$NP_DOC_POC_monthly_AVG,na.rm=TRUE)
  GPP_NPP[i,"maxGPPminusNPP"]<-max(myvals$GPPminusNPP,na.rm=TRUE)[1]
  GPP_NPP[i,"minGPPminusNPP"]<-min(myvals$GPPminusNPP,na.rm=TRUE)
  GPP_NPP[i,"meanGPPminusNPP"]<-mean(myvals$GPPminusNPP,na.rm=TRUE)}
#Export GPP/NPP parameters as excel file:
#write_xlsx(GPP_NPP, 
        #   path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/GPP_NPP.xlsx")


##Format raw data for sampling period spring 2009-winter 2010/2011
GPP_NPP.raw<-datall[which(datall$Date_Time>="2009-03-09"),]
GPP_NPP.raw<-GPP_NPP.raw[which(GPP_NPP.raw$Date_Time<="2011-02-15"),]
GPP_NPP.raw<-GPP_NPP.raw[,c("year","month","Date_Time","Biom_C","GPP_C","NP_DOC_POC")]
colnames(GPP_NPP.raw)<-c("year","month","date","CBIOMmg_per_m3","GPP","NPP")
GPP_NPP.raw$GPPminusNPP<-GPP_NPP.raw$GPP - GPP_NPP.raw$NPP
GPP_NPP.raw$NPPvsGPP<-GPP_NPP.raw$NPP / GPP_NPP.raw$GPP
GPP_NPP.raw$RESPvsGPP<-GPP_NPP.raw$GPPminusNPP /GPP_NPP.raw$GPP
GPP_NPP.raw$GPPvsBIOM<-GPP_NPP.raw$GPP / GPP_NPP.raw$CBIOMmg_per_m3 
GPP_NPP.raw$NPPvsBIOM<-GPP_NPP.raw$NPP / GPP_NPP.raw$CBIOMmg_per_m3
GPP_NPP.raw$RESPvsBIOM<-GPP_NPP.raw$GPPminusNPP / GPP_NPP.raw$CBIOMmg_per_m3
#write_xlsx(GPP_NPP.raw, 
 #  path="C:/Hannah/Biological Oceanography/Master Thesis/Project/LIM/GPP_NPP_RATIOS.xlsx")
