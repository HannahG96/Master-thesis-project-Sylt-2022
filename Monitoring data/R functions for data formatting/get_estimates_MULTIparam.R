#CALCULATION OF PARAMETER ESTIMATES-MULTIPLE PARAMETERS
#-->seasonal AVG
#-->seasonal SD & %SD
##-->annual MIN/MAX
#-->annual AVG

#FUNCTION PARAMETERS
#dat=data frame with data
#timeseries=years of the timeseries
#stations=names of sampling stations
#colnames.params=column names of parameters

get_estimates_MULTIparam<-function(dat,timeseries,stations,colnames.params){
  
#1.Monthly AVG + SD of SRB
#=monthly AVG calculated based on means of each station
#=monthly SD calculated based on pooled data across stations
  months<-1:12
  Nobs_stations<-paste("Nobs",stations,sep="_")
  colnames.df_month<-c(paste(colnames.params,"monthly_AVG",sep="_"),
                       paste(colnames.params,"monthly_SD",sep="_"))
  
  data.param_month<-as.data.frame(matrix(NA,ncol=(2+length(Nobs_stations)+length(colnames.df_month)),nrow=(length(months)*length(timeseries)),
                                         dimnames=list(NULL,c("year","month",Nobs_stations,colnames.df_month))))
  whichcols<-vector(mode="numeric",length=length(colnames.params))
  for(i in 1:length(colnames.params))whichcols[i]<-which(colnames(dat)==colnames.params[i])
  
  colstart<-length(c("year","month",Nobs_stations))+1
  nrow.df<-0
  for(i in 1:length(timeseries)){
    for(u in 1:length(months)){
      nrow.df<-nrow.df+1
      data.param_month[nrow.df,"year"]<-timeseries[i]
      data.param_month[nrow.df,"month"]<-months[u]
      mydat<-dat[which(dat[,"year"]==timeseries[i]),]
      mydat<-mydat[which(mydat[,"month"]==months[u]),]
      monthly.mean<-as.data.frame(matrix(NA,nrow=length(colnames.params),ncol=length(stations),dimnames=list(NULL,NULL)))
      for(p in 1:length(stations)){
        if(length(stations)>1){
        mydat2<-mydat[which(mydat[,"station"]==stations[p]),]
        }else{mydat2<-mydat}
        data.param_month[nrow.df, Nobs_stations[p]]<-nrow(mydat2)
        monthly.mean[,p]<-unname(apply(mydat2[,whichcols],2,FUN=mean,na.rm=TRUE))}#monthly AVG per station
  data.param_month[nrow.df,colstart:(colstart+length(colnames.params)-1)]<-
        unname(apply(monthly.mean,1,FUN=mean,na.rm=TRUE))#monthly AVG across stations
  data.param_month[nrow.df,(colstart+length(colnames.params)):(colstart+2*length(colnames.params)-1)]<-
        unname(apply(mydat[,whichcols],2,FUN=sd,na.rm=TRUE))}}#monthly SD of pooled data
  
#2.seasonal AVG + SD + %SD of SRB
#-->seasonal AVG calculated based on monthly AVG
#-->seasonal SD calculated based on pooled data across months & stations
#ORDER OF SEASONS: winter-->spring-->summer-->autumn
#START: spring of 1st year of timeseries
#END: autumn of last year of timeseries
  spring<-3:5 #March-May
  summer<-6:8 #June-August
  autumn<-9:11 #September-November
  winter<-c(12,1,2)#December(previous year)-February
  data.param_month$season<-dat$season<-NA#add season column
  for(i in 1:3){ 
    data.param_month[which(data.param_month[,"month"]==spring[i]),"season"]<-"spring"
    data.param_month[which(data.param_month[,"month"]==summer[i]),"season"]<-"summer"
    data.param_month[which(data.param_month[,"month"]==autumn[i]),"season"]<-"autumn"
    data.param_month[which(data.param_month[,"month"]==winter[i]),"season"]<-"winter"
    dat[which(dat[,"month"]==spring[i]),"season"]<-"spring"
    dat[which(dat[,"month"]==summer[i]),"season"]<-"summer"
    dat[which(dat[,"month"]==autumn[i]),"season"]<-"autumn"
    dat[which(dat[,"month"]==winter[i]),"season"]<-"winter"}
  data.param_month$yearS<-data.param_month$year#add a special year column that includes december of PREVIOUS year
  dat$yearS<-dat$year
 # for(i in 1:2){
 #   data.param_month[which(data.param_month[,"month"]==i),"yearS"]<-data.param_month[which(data.param_month[,"month"]==i),"year"]-1
 #   dat[which(dat[,"month"]==i),"yearS"]<-dat[which(dat[,"month"]==i),"year"]-1}
  data.param_month[which(data.param_month[,"month"]==12),"yearS"]<-data.param_month[which(data.param_month[,"month"]==12),"year"]+1
  dat[which(dat[,"month"]==12),"yearS"]<-dat[which(dat[,"month"]==12),"year"]+1
  
  colnames.df_season<-c(paste(colnames.params,"seasonal_AVG",sep="_"),
                       paste(colnames.params,"seasonal_SD",sep="_"),
                       paste(colnames.params,"seasonal_pctSD",sep="_"))
  data.param_season<-as.data.frame(matrix(NA,ncol=(2+length(Nobs_stations)+length(colnames.df_season)),nrow=(length(timeseries)*4),
                                          dimnames=list(NULL,c("year","season",Nobs_stations,colnames.df_season))))
  season<-c("winter","spring","summer","autumn")
  data.param_season$season<-c(rep(season,length(timeseries)))
  
  whichcols.monthly_AVG<-vector(mode="numeric",length=length(colnames.params))
  colnames.params.monthly_AVG<-paste(colnames.params,"monthly_AVG",sep="_")
  for(i in 1:length(colnames.params))whichcols.monthly_AVG[i]<-which(colnames(data.param_month)==colnames.params.monthly_AVG[i])
  
  colstart<-length(c("year","season",Nobs_stations))+1
  
  nrow.df<-0
  for(i in 1:length(timeseries)){
    #if(i<length(timeseries)){
      x<-3
    #}else{ #if(i==length(timeseries)){
     # x<-2
      #season<-c("spring","summer","autumn")}
    for(u in 1:length(season)){
      nrow.df<-nrow.df+1
      data.param_season[nrow.df,"year"]<-timeseries[i]
      mydat<-data.param_month[which(data.param_month[,"yearS"]==timeseries[i]),]
      mydat<-mydat[which(mydat[,"season"]==season[u]),]
      for(p in 1:length(stations)){
      data.param_season[nrow.df,Nobs_stations[p]]<-sum(mydat[,Nobs_stations[p]])}
      data.param_season[nrow.df,colstart:(colstart+length(colnames.params)-1)]<-
        unname(apply(mydat[,whichcols.monthly_AVG],2,FUN=mean,na.rm=TRUE))#seasonal AVG calculated based on monthly AVG
      mydat2<-dat[which(dat[,"yearS"]==timeseries[i]),]
      mydat2<-mydat2[which(mydat2[,"season"]==season[u]),]
      data.param_season[nrow.df,(colstart+length(colnames.params)):(colstart+2*length(colnames.params)-1)]<-
        unname(apply(mydat2[,whichcols],2,FUN=sd,na.rm=TRUE))
      data.param_season[nrow.df,(colstart+2*length(colnames.params)):(colstart+3*length(colnames.params)-1)]<-
        as.numeric(data.param_season[nrow.df,(colstart+length(colnames.params)):(colstart+2*length(colnames.params)-1)])*100/as.numeric(data.param_season[nrow.df,colstart:(colstart+length(colnames.params)-1)])
      }}
  
#3.Annual AVG + MIN + MAX of SRB
#--> annual AVG calculated based monthly AVG
#-->annual MIN/MAX calculated based on pooled data
  Nobs.sum_stations<-paste("Nobs.sum",stations,sep="_")
  colnames.df_annual<-c(paste(colnames.params,"annual_AVG",sep="_"),
                        paste(colnames.params,"annual_MAX",sep="_"),
                        paste(colnames.params,"annual_MIN",sep="_"))
  data.param_annual<-as.data.frame(matrix(NA,ncol=(1+length(Nobs.sum_stations)+length(colnames.df_annual)),nrow=length(timeseries),
                                          dimnames=list(NULL,c("year",Nobs.sum_stations,colnames.df_annual))))
  colstart<-length(c("year",Nobs.sum_stations))+1
  nrow.df<-0
  for(i in 1:length(timeseries)){
    nrow.df<-nrow.df+1
    data.param_annual[nrow.df,"year"]<-timeseries[i]
    mydat<-data.param_month[which(data.param_month[,"year"]==timeseries[i]),]
    data.param_annual[nrow.df,colstart:(colstart+length(colnames.params)-1)]<-
      unname(apply(mydat[,whichcols.monthly_AVG],2,FUN=mean,na.rm=TRUE))#annual AVG
    for(u in 1:length(stations))data.param_annual[nrow.df,Nobs.sum_stations[u]]<-sum(mydat[,Nobs_stations[u]])
    mydat2<-dat[which(dat[,"year"]==timeseries[i]),]
    data.param_annual[nrow.df,(colstart+length(colnames.params)):(colstart+2*length(colnames.params)-1)]<-
      unname(apply(mydat2[,whichcols],2,FUN=max,na.rm=TRUE))#annual MAX
    data.param_annual[nrow.df,(colstart+2*length(colnames.params)):(colstart+3*length(colnames.params)-1)]<-
      unname(apply(mydat2[,whichcols],2,FUN=min,na.rm=TRUE))}#annual MIN
  
#4.Restructure monthly/seasonal/annual data frames for function output 
  minus<-c(which(colnames(data.param_month)=="season"),which(colnames(data.param_month)=="yearS"))
  data.param_month<-data.param_month[,-minus]
  data.param<-merge(data.param_annual,data.param_season,by="year",all=TRUE)
  output<-list(data.param_month,data.param)
  names(output)<-c("Monthly estimates","Seasonal + Annual estimates")
  
  return(output)
}#END OF FUNCTION
