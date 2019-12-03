rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(33)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ncomm=10
ngibbs=1000
nburn=ngibbs/2

#priors
psi=0.01
gamma=0.1
#----------------------------------------------------------
#run LDA no covariates to get initial values

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
source('LDA.abundance main function.R')
sourceCpp('aux1.cpp')

res=LDA.abundance(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,psi=psi,gamma=gamma)

nloc=nrow(y)
nspp=ncol(y)
array.lsk.init=res$array.lsk

#determine optimal number of groups
nlk=apply(array.lsk.init,c(1,3),sum)
theta=nlk/apply(nlk,1,sum)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(theta)
ncomm=9

prop=apply(theta>0.8,2,sum,na.rm=T) #see which communities are never above 0.8
which(prop!=0)
cond=prop!=0
ncomm=sum(cond)

#re-distribute individuals within array.lsk.init that are in eliminated communities
array.lsk=array.lsk.init[,,cond]
for (i in 1:nloc){
  for (j in 1:nspp){
    tmp=array.lsk.init[i,j,!cond]
    n=sum(tmp)
    if (n>0){
      z=rmultinom(1,size=n,prob=rep(1/ncomm,ncomm))
      array.lsk[i,j,]=array.lsk[i,j,]+z
    }
  }
}

#export results
dat1=matrix(array.lsk,nloc*nspp*ncomm,1)
setwd('U:\\GIT_models\\git_LDA_MS')
write.csv(dat1,'array lsk.csv',row.names=F)