# rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs sampler main function.R')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data8.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat8.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ncomm=10
ngibbs=1000
nburn=ngibbs/2
phi.prior=0.1
a1=b1=1
b2=20
a2=10

res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior,a1=a1,b1=b1,a2=a2,b2=b2)

plot(res$lambda[ngibbs,],type='h')
plot(res$llk,type='l')