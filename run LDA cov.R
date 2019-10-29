# rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(10)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ncomm=5
phi.prior=0.01
a.gamma=b.gamma=0.1

#----------------------------------------
#run LDA no covariates to get initial values

ngibbs=1000
nburn=ngibbs/2

#get functions
setwd('U:\\GIT_models\\LdaPoisson_nocov')
source('LdaPoisson_nocov main function.R')
source('LdaPoisson_nocov aux functions.R')
sourceCpp('LdaPoisson_nocov_aux_cpp.cpp')

res=gibbs.LDA.nocov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,
                    phi.prior=phi.prior,a.gamma=a.gamma,b.gamma=b.gamma)

nloc=nrow(y)
nspp=ncol(y)
array.lsk.init=array(res$array.lsk[ngibbs,],dim=c(nloc,nspp,ncomm))
#----------------------------------------
#LDA with covariates

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA cov main function.R')
source('LDA cov aux functions.R')
sourceCpp('LDA_cov_aux1_cpp.cpp')

var.betas=10
res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                  phi.prior=phi.prior,array.lsk.init=array.lsk.init,
                  a.gamma=a.gamma,b.gamma=b.gamma,var.betas=var.betas)

plot(res$lambda[ngibbs,],type='h')
plot(res$llk,type='l')