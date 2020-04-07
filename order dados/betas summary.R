rm(list=ls(all=TRUE))

#get data
setwd('Z:\\Users\\drvalle\\GIT_models\\git_LDA_MS\\order dados')
dat=read.csv('y.csv',as.is=T)
xmat=data.matrix(read.csv('xmat.csv',as.is=T))
nparam=ncol(xmat)
colnames(xmat)=c('intercept','greater6days','saturday','tueswed')

#get parameters
betas=read.csv('betas.csv',as.is=T)
ncomm=ncol(betas)/nparam; ncomm
tmp=apply(betas,2,quantile,c(0.025,0.5,0.975))
betas.med=matrix(tmp[2,],nparam,ncomm)
rownames(betas.med)=colnames(xmat)
grupos=paste0('gr',1:ncomm)
colnames(betas.med)=grupos

#get 95% CI
betas.lo=matrix(tmp[1,],nparam,ncomm)
betas.hi=matrix(tmp[3,],nparam,ncomm)
cond.mat=betas.lo<0 & betas.hi>0
signif1=!cond.mat #statistically significant results
rownames(signif1)=colnames(xmat)
colnames(signif1)=grupos

#output results
write.csv(signif1,'betas signif.csv')
write.csv(betas.med,'betas summary.csv')