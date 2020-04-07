rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
set.seed(10)

#get data
setwd('Z:\\Users\\drvalle\\GIT_models\\git_LDA_MS\\order dados')
dat=read.csv('y.csv',as.is=T)
xmat=data.matrix(read.csv('xmat.csv',as.is=T))
y=data.matrix(dat)
nloc=nrow(y)
nspp=ncol(y)

#get array.lsk
tmp=read.csv('array lsk.csv',as.is=T)
ncomm=length(tmp$V1)/(nspp*nloc); ncomm
array.lsk.init=array(tmp$V1,dim=c(nloc,nspp,ncomm))

#basic settings
ngibbs=10000
nburn=ngibbs/2

#priors
phi.prior=0.01
var.betas=c(10,rep(10,ncol(xmat)-1))
gamma=0.1

#----------------------------------------------------------
#LDA with covariates

#get functions
setwd('Z:\\Users\\drvalle\\GIT_models\\git_LDA_MS')
source('LDA cov main function.R')
source('LDA cov aux functions.R')
sourceCpp('LDA_cov_aux1_cpp.cpp')
sourceCpp('slice_betas.cpp')
sourceCpp('slice_NBN.cpp')

res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior,array.lsk.init=array.lsk.init,
                  var.betas=var.betas)
plot(res$llk,type='l')
plot(res$fmodel,type='l')

seq1=8500:ngibbs
par(mfrow=c(4,3),mar=rep(1,4))
for (i in 1:12){
  plot(res$betas[seq1,i],type='l')
}

setwd('Z:\\Users\\drvalle\\GIT_models\\git_LDA_MS\\order dados')
write.csv(res$phi[seq1,],'phi.csv',row.names=F)
write.csv(res$nlk[seq1,],'nlk.csv',row.names=F)
write.csv(res$betas[seq1,],'betas.csv',row.names=F)
write.csv(res$NBN[seq1],'NBN.csv',row.names=F)