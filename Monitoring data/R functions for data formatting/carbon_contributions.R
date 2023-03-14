##################################################################
#MEAN CARBON CONTRIBUTION (based on average of contribution of each species/group to respective carbon pool at each sampling day)
#meanCC=mean(carbon of species i at sampling day x/carbon pool at sampling day x)
#-->See method section about "Weighted Community Biomass" of Julien Meunier
#-->based on seasonal average daily contribution of species to compartment carbon
##################################################################

source("R functions/get_estimates_ONEparam.R")
#parameters:
#dat=dataframe with species carbons
#stations=name(s) of sampling stations
#species=list of species names per compartment (each list element represent on compartment)
#columns=list of column names in dat of species per compartment (each list element represent on compartment)

spptocompC<-function(dat,stations,species,columns,comps){

#Extract df with total microplankton carbon of SRB: Dec. 2008-Feb.2011
TOTCARB<-dat[which(dat[,"date"]>="2008-12-01"),]
TOTCARB<-TOTCARB[which(TOTCARB[,"date"]<="2011-02-28"),]
#Create season column:
TOTCARB$season<-paste(TOTCARB$year,TOTCARB$month,sep=" ")
addDec1<-which(TOTCARB$season=="2008 12")
addDec2<-which(TOTCARB$season=="2009 12")
addDec3<-which(TOTCARB$season=="2010 12")
for(i in 1:nrow(TOTCARB)){
  if(length(which(addDec1==i))>0||length(which(addDec2==i))>0||length(which(addDec3==i))>0){
    TOTCARB[i,"season"]<-paste((TOTCARB[i,"year"]+1),"winter",sep=" ")  
  }else if(TOTCARB[i,"month"]==1 || TOTCARB[i,"month"]==2){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"winter",sep=" ")
  }else if(TOTCARB[i,"month"]==3 || TOTCARB[i,"month"]==4 || TOTCARB[i,"month"]==5){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"spring",sep=" ")
  }else if(TOTCARB[i,"month"]==6 || TOTCARB[i,"month"]==7 || TOTCARB[i,"month"]==8){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"summer",sep=" ")
  }else if(TOTCARB[i,"month"]==9 || TOTCARB[i,"month"]==10 || TOTCARB[i,"month"]==11){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"autumn",sep=" ")}}

#Create df that indicates carbon contributions
CARBcontr<-as.data.frame(matrix(NA,nrow=0,ncol=4,
                                dimnames=list(NULL,
                                              c("Species","Compartment","Season","pct_Carbon"))))

#####CARBON CONTRIBUTION OF EACH INCLUDED SPP TO COMPARTMENT CARBON
##############ACROSS SEASONS 2009/2010 (including Dec 2008 & Jan.+Feb. 2011)
n<-1
for(u in 1:length(comps)){
#Total compartment carbon:
TOTCARB$Comp<-apply(as.data.frame(TOTCARB[,columns[[u]]]),1,sum,na.rm=TRUE)
#species carbon to compartment carbon:
for(i in 1:length(species[[u]])){
if(length(stations)==1){
  df<-data.frame(year=TOTCARB$year,
                 month=TOTCARB$month,
                 Date=TOTCARB$date,
                 CarbonContr=TOTCARB[,which(colnames(TOTCARB)==columns[[u]][i])]/TOTCARB$Comp)
}else{
  df<-data.frame(year=TOTCARB$year,
                 month=TOTCARB$month,
                 Date=TOTCARB$date,
                 CarbonContr=TOTCARB[,which(colnames(TOTCARB)==columns[[u]][i])]/TOTCARB$Comp,
                 station=TOTCARB$station)}
  df<-get_estimates_ONEparam(dat=df,timeseries=2008:2011,
                             stations=stations,colname.param="CarbonContr")[[2]]
  df$season<-paste(df$year,df$season,sep=" ")
  
  myspp<-species[[u]][i]
  for(p in 1:length(unique(TOTCARB$season))){
    CARBcontr[n,"Species"]<-myspp
    CARBcontr[n,"Compartment"]<-comps[u]
    CARBcontr[n,"Season"]<-unique(TOTCARB$season)[p]
    CARBcontr[n,"pct_Carbon"]<-
      df[which(df[,"season"]==unique(TOTCARB$season)[p]),"seasonal_AVG"]*100
    n<-n+1}}}

return(CARBcontr)
}#END OF FUNCTION

##########################################################################
#B)COMPARTMENT CARBON DISTRIBUTION TO TOTAL SRB CARBON POOL OF FUNCTIONAL GROUP

#parameters:
#dat=dataframe with compartment carbons & total SRB functional group carbons
#stations=name(s) of sampling stations
#comps=vector of compartment names 
#cols.comps=vector of column names of respective compartments
#totC=vector of functional group names representing each SRB carbon pool
#cols.totC=vector of column names of respective SRB carbon pool

comptototC<-function(dat,stations,comps,cols.comps,totC,cols.totC){

#Extract df with total microplankton carbon of SRB: Dec. 2008-Feb.2011
TOTCARB<-dat[which(dat[,"date"]>="2008-12-01"),]
TOTCARB<-TOTCARB[which(TOTCARB[,"date"]<="2011-02-28"),]
#Create season column:
TOTCARB$season<-paste(TOTCARB$year,TOTCARB$month,sep=" ")
addDec1<-which(TOTCARB$season=="2008 12")
addDec2<-which(TOTCARB$season=="2009 12")
addDec3<-which(TOTCARB$season=="2010 12")
for(i in 1:nrow(TOTCARB)){
  if(length(which(addDec1==i))>0||length(which(addDec2==i))>0||length(which(addDec3==i))>0){
    TOTCARB[i,"season"]<-paste((TOTCARB[i,"year"]+1),"winter",sep=" ")  
  }else if(TOTCARB[i,"month"]==1 || TOTCARB[i,"month"]==2){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"winter",sep=" ")
  }else if(TOTCARB[i,"month"]==3 || TOTCARB[i,"month"]==4 || TOTCARB[i,"month"]==5){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"spring",sep=" ")
  }else if(TOTCARB[i,"month"]==6 || TOTCARB[i,"month"]==7 || TOTCARB[i,"month"]==8){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"summer",sep=" ")
  }else if(TOTCARB[i,"month"]==9 || TOTCARB[i,"month"]==10 || TOTCARB[i,"month"]==11){
    TOTCARB[i,"season"]<-paste(TOTCARB[i,"year"],"autumn",sep=" ")}}

#Create df that indicates carbon contributions
CARBcontr<-as.data.frame(matrix(NA,nrow=0,ncol=4,
                                dimnames=list(NULL,
                c("Compartment","CarbonPool","Season","pct_Carbon"))))

#####CARBON CONTRIBUTION OF EACH COMPARTMENT TO TOTAL FUNCTIONAL GROUP CARBON
##############ACROSS SEASONS 2009/2010 (including Dec 2008 & Jan.+Feb. 2011)
n<-1
for(u in 1:length(comps)){
if(length(stations)==1){
    df<-data.frame(year=TOTCARB$year,
                   month=TOTCARB$month,
                   Date=TOTCARB$date,
                   CarbonContr=TOTCARB[,which(colnames(TOTCARB)==cols.comps[u])]/
                   TOTCARB[,which(colnames(TOTCARB)==cols.totC[u])])
}else{
  df<-data.frame(year=TOTCARB$year,
                 month=TOTCARB$month,
                 Date=TOTCARB$date,
                 CarbonContr=TOTCARB[,which(colnames(TOTCARB)==cols.comps[u])]/
                 TOTCARB[,which(colnames(TOTCARB)==cols.totC[u])],
                 station=TOTCARB$station)}
    df<-get_estimates_ONEparam(dat=df,timeseries=2008:2011,
                               stations=stations,colname.param="CarbonContr")[[2]]
    df$season<-paste(df$year,df$season,sep=" ")

    for(p in 1:length(unique(TOTCARB$season))){
      CARBcontr[n,"Compartment"]<-comps[u]
      CARBcontr[n,"CarbonPool"]<-totC[u]
      CARBcontr[n,"Season"]<-unique(TOTCARB$season)[p]
      CARBcontr[n,"pct_Carbon"]<-
        df[which(df[,"season"]==unique(TOTCARB$season)[p]),"seasonal_AVG"]*100
      n<-n+1}}

return(CARBcontr)
}#END OF FUNCTION

#########################################################################
####PLOT OF RESULTS

plot.Ccontributions<-function(CARBcontr,givetitle=NA){
  
if(is.na(givetitle)==TRUE)givetitle<-""
colnames(CARBcontr)[1:2]<-c("whatfill","whatfacet")
CARBcontr$Season<-factor(CARBcontr$Season, 
     levels=c("2009 winter","2009 spring","2009 summer","2009 autumn",
              "2010 winter","2010 spring","2010 summer","2010 autumn",
              "2011 winter"))
CARBcontr$whatfill<-factor(CARBcontr$whatfill,levels=unique(CARBcontr$whatfill))
CARBcontr$whatfacet<-factor(CARBcontr$whatfacet,levels=unique(CARBcontr$whatfacet))
cols<-c("gold","turquoise3","forestgreen","royalblue","mediumorchid","red","yellow","pink",
         "lawngreen","chocolate3","magenta","limegreen","aquamarine4","darkgoldenrod4",
         "firebrick3","lightseagreen","mediumvioletred","chartreuse","darkred","darkslateblue",
        "lightcoral","khaki","sandybrown","lightsteelblue","palegreen","thistle",
        "yellow3","tan1","plum","grey","salmon","lightskyblue3","olivedrab1","orange",
        "mediumaquamarine","seashell3","rosybrown","burlywood3","darkolivegreen4","cornflowerblue")
if(length(unique(CARBcontr$whatfill))<=length(cols)){
cols.plot<-sample(cols,length(unique(CARBcontr$whatfill)),replace=FALSE)
}else{
cols.plot<-c(sample(cols,length(cols),replace=FALSE),
       rep("gray21",(length(unique(CARBcontr$whatfill))-length(cols))))}

plotC<-
ggplot(CARBcontr, aes(fill=whatfill, y=pct_Carbon, x=Season))+ 
  theme_minimal()+
  geom_bar(position = "stack", stat="identity", width = 0.5, color="#90A4ADFF")+
  scale_fill_manual("",values=cols.plot,
                    breaks=unique(CARBcontr$whatfill),
                    labels=unique(CARBcontr$whatfill))+
  ggtitle(givetitle)+
  ylab("% of Total Carbon")+
  theme(axis.text.x=element_text(angle=20,size=9,face="bold",color="navy"),
        axis.text.y=element_text(face="bold",color="slategray3"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(face="bold", size=8,color="slategray3"),
        panel.grid.major = element_line(colour="lightgrey", linetype="dashed"), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.2, 'cm'))+
    facet_wrap(~whatfacet,ncol=1)
return(plotC)
}#END OF FUNCTION
