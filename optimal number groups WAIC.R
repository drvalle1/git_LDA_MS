rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(10)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA cov and LDA no cov.R')
dat=read.csv('fake data.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ngibbs=1000
nburn=ngibbs/2

#priors
phi.prior=0.01
a.gamma=b.gamma=0.1
var.betas=10

waic=rep(NA,15)
for (ncomm in 3:15){
  #run Gibbs sampler
  res=LDAcov(y=y,xmat=xmat,ncomm=ncomm,
             phi.prior=phi.prior,a.gamma=a.gamma,b.gamma=b.gamma,
             ngibbs=ngibbs,nburn=nburn,var.betas=var.betas)
  
  #calculate WAIC
  # plot(res$llk,type='l')
  seq1=nburn:ngibbs
  tmp=res$llk.ind.out[seq1,]
  p1=-2*sum(log(colMeans(exp(tmp))))
  p2=2*sum(apply(tmp,2,var))
  waic[ncomm]=p1+p2;
}

par(mfrow=c(1,1))
plot(waic,type='l')