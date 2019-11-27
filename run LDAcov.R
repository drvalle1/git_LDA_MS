rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(10)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
y=data.matrix(dat)
nloc=nrow(y)
nspp=ncol(y)

#get array.lsk
tmp=read.csv('array lsk.csv',as.is=T)
ncomm=5 #this was estimated based on a previous step
array.lsk=array(tmp$V1,dim=c(nloc,nspp,ncomm))

#basic settings
ngibbs=1000
nburn=ngibbs/2

#priors
phi.prior=0.01
a.gamma=b.gamma=0.1
var.betas=10
gamma=0.1

#----------------------------------------------------------
#LDA with covariates

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA cov main function.R')
source('LDA cov aux functions.R')
sourceCpp('LDA_cov_aux1_cpp.cpp')

res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior,array.lsk.init=array.lsk,
                  a.gamma=a.gamma,b.gamma=b.gamma,var.betas=var.betas)
plot(res$llk,type='l')
