#EXTRACT SPECIES NAMES LIST OF SAMPLED SPECIMENS FROM A TAB FILE
get_species<-function(myfile){
read.lines<-0
b<-"\tDEPTH"
repeat{
  current.line<-read.lines+1
  a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
  aa<-strsplit(a,split=",",fixed=TRUE)[[1]][1]
  read.lines<-current.line
  if(aa==b){break}}
spp.names<-col.names<-vector()
repeat{
  current.line<-read.lines+1
  a<-scan(file=myfile,sep="\n",skip=read.lines,nlines=1,what="character")
  a<-sub(", biomass as carbon","",a,fixed=TRUE)
  aa<-strsplit(a,split="[",fixed=TRUE)
  if(length(aa[[1]])==1){break
  }else{
  spp.names<-c(spp.names,aa[[1]][1])
  b<-strsplit(a,split=")",fixed=TRUE)[[1]][1]
  bb<-strsplit(b,split="(",fixed=TRUE)[[1]][2]
  col.names<-c(col.names,bb)
  read.lines<-current.line}}
spp.names<-sub("\t","",spp.names,fixed=TRUE)
colnames.list<-strsplit(col.names,split=" ", fixed=TRUE)
column.names<-vector()
for(i in 1:length(colnames.list)){
  mycol<-colnames.list[[i]]
  mycol<-sub("/","_",x=mycol,fixed=TRUE)
  mycol<-sub("+","_",x=mycol,fixed=TRUE)
  for(u in 1:length(mycol)){
    if(u==1){mycolname<-mycol[u]
    }else{mycolname<-paste(mycolname,mycol[u],sep="_")}}
  column.names<-c(column.names,mycolname)}

return(data.frame(species=spp.names,column=column.names))
}
