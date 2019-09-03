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
ncomm=8
ngibbs=1000
nburn=ngibbs/2
phi.prior=0.1
lambda.a=lambda.b=0.1

res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior,lambda.a=lambda.a,lambda.b=lambda.b)


