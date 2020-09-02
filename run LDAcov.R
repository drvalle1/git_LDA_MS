rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
set.seed(10)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
# xmat=xmat[,1:2]
y=data.matrix(dat)
nloc=nrow(y)
nspp=ncol(y)

#get array.lsk
tmp=read.csv('array lsk.csv',as.is=T)
ncomm.init=length(tmp$V1)/(nspp*nloc); ncomm.init
array.lsk.init=array(tmp$V1,dim=c(nloc,nspp,ncomm.init))

#get phi
phi.init=data.matrix(read.csv('phi step1.csv',as.is=T))

#basic settings
ngibbs=1000
nburn=ngibbs/2

#priors
phi.prior=0.01
var.betas=c(10,rep(10,ncol(xmat)-1))
gamma=0.1

#----------------------------------------------------------
#determine optimal number of groups

nlk=apply(array.lsk.init,c(1,3),sum)
theta1=nlk/apply(nlk,1,sum)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(theta1)

#----------------------------------------------------------
#LDA with covariates

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA cov main function.R')
source('LDA cov aux functions.R')
sourceCpp('LDA_cov_aux1_cpp.cpp')
sourceCpp('slice_betas.cpp')
sourceCpp('slice_NBN.cpp')

#re-estimate phi
ncomm=4
res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,y=y,xmat=xmat,
                  phi.prior=phi.prior,array.lsk.init=array.lsk.init,
                  var.betas=var.betas,estimate.phi=T,
                  phi.init=phi.init)

#do not re-estimate phi
res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,y=y,xmat=xmat,
                  phi.prior=phi.prior,array.lsk.init=array.lsk.init,
                  var.betas=var.betas,phi.init=phi.init,estimate.phi=F)