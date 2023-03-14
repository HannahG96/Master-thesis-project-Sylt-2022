######FUNCTION TO IMPORT PANGAEA TAB FILES
#-->imports files stored in "C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly"
#-->input "myfile" should also contain the folder name where data is stored within "Data assembly" e.g. "myfolder/myfile.tab)
#-->all columns of the output data frame are of class "character" (columns need to be manually converted into "numeric"), if needed)

import_pangaea_file<-function(myfile){

myfile<-paste("C:/Hannah/Biological Oceanography/Master Thesis/Project/Data assembly/",myfile,sep="")

###Determine start and end line of metadata
read.lines<-0
b<-"*/"
repeat{
  current.line<-read.lines+1
  a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
  read.lines<-current.line
  if(a==b){break}}
metadata.end<-read.lines
#-->data starts at read.lines+1

###Determine number and names of columns
#Start of parameter list:
read.lines<-0
b<-"Parameter(s):"
repeat{
  current.line<-read.lines+1
  a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
  aa<-strsplit(a,split="\t",fixed=TRUE)[[1]][1]
  read.lines<-current.line
  if(aa==b){break}}
param.start<-read.lines
#-->list of parameters starts at read.lines

#Extract names of all parameters:
read.lines<-param.start-1
repeat{
  if(read.lines==(param.start-1)){
a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
aa<-strsplit(a,split=")",fixed=TRUE)[[1]][2]
aaa<-strsplit(aa,split="(",fixed=TRUE)[[1]][2]
param.names<-aaa
}else{
a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
if(strsplit(a,split=":",fixed=TRUE)[[1]][1]=="License")break
aa<-strsplit(a,split=")",fixed=TRUE)[[1]][1]
aaa<-strsplit(aa,split="(",fixed=TRUE)[[1]][2]
param.names<-c(param.names,aaa)}

read.lines<-read.lines+1
}

###Import data table
#Determine start of data table (=numbers)
read.lines<-metadata.end
params<-c()
repeat{
if(length(params)==length(param.names))break
a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
aa<-strsplit(a,split=c("\t"),fixed=TRUE)[[1]]
params<-c(params,aa)
read.lines<-read.lines+1}
data.start<-read.lines+1 #-->data table starts at read.lines+1

#Rename parameters for column names of data table:
colnames.list<-strsplit(param.names,split=" ", fixed=TRUE)
colnames<-vector()
for(i in 1:length(colnames.list)){
mycol<-colnames.list[[i]]
mycol<-sub("/","_",x=mycol,fixed=TRUE)
mycol<-sub("+","_",x=mycol,fixed=TRUE)
#mycol<-sub("sp.","sp",x=mycol,fixed=TRUE)
for(u in 1:length(mycol)){
  if(u==1){mycolname<-mycol[u]
  }else{mycolname<-paste(mycolname,mycol[u],sep="_")}}
colnames<-c(colnames,mycolname)}

#Import all rows of data table and store values in data frame:
dat<-as.data.frame(matrix(NA,nrow=0, ncol=length(param.names)))#df to store data values
repeat{
a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
if(length(a)==0)break
aa<-strsplit(a,split="\t",fixed=TRUE)#character vector
dat<-rbind(dat,aa[[1]])
read.lines<-read.lines+1}
if(ncol(dat)==length(colnames)){
colnames(dat)<-colnames
}else{
  print("Error in reading data table: colum number unequal to parameter number")
}

#return the data frame:
print(paste("Identify",length(colnames),"parameters",sep=" "))
print(paste("Start reading-in data at row",data.start, "of tab file",sep=" "))
print(paste("Read",read.lines,"lines of tab file",sep=" "))
return(dat)#Note: numeric values of data table are in character format

} #END OF FUNCTION


