rm(list=ls(all=TRUE))
library('Rcpp')
library('mvtnorm')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
source('gibbs sampler main function.R')
sourceCpp('aux1.cpp')

#get data
tmp=read.csv('fake data4.csv',as.is=T)
ind=which(colnames(tmp)=='X')
dat=tmp[,-ind]
nloc=nrow(dat)

dmat=read.csv('fake data4 dmat.csv',as.is=T)
dmat=data.matrix(dmat)

#basic settings
mu=3
ncomm=4
ngibbs=1000
phi.prior=0.1
nburn=ngibbs/2

#fit model
res=lda.abundance.regression(dat=dat,ncomm=ncomm,phi.prior=phi.prior,ngibbs=ngibbs,mu=mu,nburn=nburn,dmat=dmat)
