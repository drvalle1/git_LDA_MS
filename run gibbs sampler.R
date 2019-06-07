rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
source('gibbs sampler main function.R')
sourceCpp('aux1.cpp')

#get data
tmp=read.csv('fake data8.csv',as.is=T)
ind=which(colnames(tmp)=='X')
dat=tmp[,-ind]
nloc=nrow(dat)

#basic settings
gamma1=0.1
ncomm=10
ngibbs=1000
phi.prior=0.1
nburn=ngibbs/2
mu0=1.33#2
sig2=1 #1, sqrt(0.1) (very bad: finds 4 groups), sqrt(1) (almost good: finds 6 groups),sqrt(3.72)
sd0=sqrt(sig2/nloc*4)

#fit model
res=lda.abundance.regression(dat=dat,ncomm=ncomm,phi.prior=phi.prior,gamma1=gamma1,ngibbs=ngibbs,sd0=sd0,mu0=mu0,nburn=nburn,sig2=sig2)
