rm(list=ls(all=TRUE))
library('mvtnorm')
library('Rcpp')
library('MCMCpack')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA MS functions.R')
source('LDA MS gibbs sampler.R')
sourceCpp('LDA_MS_c.cpp')

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data6.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)

#get covariates
xmat=data.matrix(read.csv('fake data cov6.csv',as.is=T))

#run gibbs sampler
ncomm=10
ngibbs=1000
nburnin=ngibbs/2
res=LDA.MS.gibbs(y=y,xmat=xmat,ncomm=ncomm,ngibbs=ngibbs,nburnin=nburnin)