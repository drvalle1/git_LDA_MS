rm(list=ls(all=TRUE))
library('mvtnorm')
library('Rcpp')
library('MCMCpack')
set.seed(4)

#-----------------------------------------------------------
#*****Unstructured LDA*****

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
source('LDA.abundance main function.R')
sourceCpp('aux1.cpp')

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data6.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)

ncomm=10
gamma=0.1
ngibbs=1000
nburn=ngibbs/2
psi=0.01
res=LDA.abundance(y=y,ncomm=ncomm,gamma=gamma,
                  ngibbs=ngibbs,nburn=nburn,psi=psi)

#get MAP estimates 
ind=which(res$llk==max(res$llk))
vmat=matrix(res$vmat[ind,],nrow(y),ncomm)
vmat1=vmat[,-ncomm]
theta=matrix(res$theta[ind,],nrow(y),ncomm) 
phi=matrix(res$phi[ind,],ncomm,ncol(y))

#are the thetas similar to what they should be?
# boxplot(theta)
# theta1=theta[,1:6]
# res=matrix(NA,6,6)
# for (i in 1:6){
#   for (j in 1:6){
#     res[i,j]=cor(cbind(theta.true[,i],theta1[,j]))[1,2]
#   }
# }
# ind=numeric()
# for (i in 1:6) ind=c(ind,which(res[i,]==max(res[i,])))
# plot(theta1[,ind],theta.true)
#-----------------------------------------------------------
#*****LDA covariates*****

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA MS functions.R')
source('LDA MS gibbs sampler.R')
sourceCpp('LDA_MS_c.cpp')

#get covariates
xmat=data.matrix(read.csv('fake data cov6.csv',as.is=T))

#run gibbs sampler
ncomm=10
ngibbs=1000
nburnin=ngibbs/2
res=LDA.MS.gibbs(y=y,xmat=xmat,ncomm=ncomm,ngibbs=ngibbs,nburnin=nburnin,
                 vmat.init=vmat1,phi.init=phi)