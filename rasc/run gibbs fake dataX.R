rm(list=ls(all=TRUE))
library('mvtnorm')
library('Rcpp')
library('MCMCpack')
set.seed(4)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
xmat.centered=data.matrix(read.csv('fake data cov.csv',as.is=T))

#----------------------------------------------
#run standard LDA

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
source('LDA.abundance main function.R')
sourceCpp('aux1.cpp')

ncomm=5
ngibbs=1000
nburn=ngibbs/2
psi=0.01
gamma=0.1
res=LDA.abundance(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,psi=psi,gamma=gamma)
phi.init=matrix(res$phi[nrow(res$phi),],ncomm,ncol(y))
theta.init=matrix(res$theta[nrow(res$theta),],nrow(y),ncomm)
vmat.init=matrix(res$vmat[nrow(res$vmat),],nrow(y),ncomm)
boxplot(theta.init) 

#Don't restrict the maximum number of groups according to theta.init because these are very different models.
#As a result, the "structured" LDA might have the need for a greater number of groups than what is suggested by the "unstructured" LDA
#-------------------------------------------------
#run "structured" LDA

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA MS functions.R')
source('LDA MS gibbs sampler.R')
sourceCpp('LDA_MS_c.cpp')

#run gibbs sampler
ngibbs=1000
nburnin=ngibbs/2
phi.prior=psi

res=LDA.MS.gibbs(y=y,xmat=xmat.centered,ncomm=ncomm,ngibbs=ngibbs,nburnin=nburnin,phi.prior=phi.prior,
                 vmat.init1=vmat.init,phi.init=phi.init,
                 mu.betas=rep(0,ncol(xmat.centered)),
                 var.betas=rep(10,ncol(xmat.centered)))


