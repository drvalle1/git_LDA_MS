rm(list=ls(all=TRUE))

#get data
setwd('Z:\\Users\\drvalle\\GIT_models\\git_LDA_MS\\order dados')
dat=read.csv('y.csv',as.is=T)
nspp=ncol(dat)

#get phi
phi=read.csv('phi.csv',as.is=T)
ncomm=ncol(phi)/nspp; ncomm
phi1=matrix(colMeans(phi),ncomm,nspp)
colnames(phi1)=colnames(dat)

#get signature items
res=numeric()
for (i in 1:nspp){
  max1=max(phi1[,i])
  ind=which(phi1[,i]==max1)
  max2=max(phi1[-ind,i])
  if (max1/max2 > 2){
    tmp=c(colnames(phi1)[i],ind)
    res=rbind(res,tmp)
  }
}
colnames(res)=c('item','grp')

#organize results
res1=data.frame(res,stringsAsFactors = F)
res1$grp=as.numeric(res1$grp)
res2=res1[order(res1$grp),]

#export
write.csv(res2,'phi summary.csv',row.names=F)