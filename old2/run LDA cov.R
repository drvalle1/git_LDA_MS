# rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(10)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA cov main function.R')

#get data
dat=read.csv('fake data8.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat8.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ncomm=5
ngibbs=1001
nburn=ngibbs/2
phi.prior=0.01

res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior)

plot(res$lambda[ngibbs,],type='h')
plot(res$llk,type='l')