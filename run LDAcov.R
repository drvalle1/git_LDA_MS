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
ncomm=10
ngibbs=1000
nburn=ngibbs/2

#priors
phi.prior=0.01
a.gamma=b.gamma=0.1
var.betas=10

#run Gibbs sampler
res=LDAcov(y=y,xmat=xmat,ncomm=ncomm,
           phi.prior=phi.prior,a.gamma=a.gamma,b.gamma=b.gamma,
           ngibbs=ngibbs,nburn=nburn,var.betas=var.betas)

#calculate WAIC
plot(res$llk,type='l')
seq1=(ngibbs/2):ngibbs
tmp=res$llk.ind.out[seq1,]
p1=-2*sum(log(colMeans(exp(tmp))))
p2=2*sum(apply(tmp,2,var))
waic=p1+p2; waic

