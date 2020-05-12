rm(list=ls(all=TRUE))
library(MCMCpack)
library('Rcpp')
library('RcppArmadillo')
set.seed(33)

#get data
setwd('U:\\GIT_models\\git_LDA_MS')
dat=read.csv('fake data.csv',as.is=T)
y=data.matrix(dat)

#basic settings
ncomm.init=10
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

res=LDA.abundance(y=y,ncomm=ncomm.init,ngibbs=ngibbs,nburn=nburn,psi=psi,gamma=gamma)

nloc=nrow(y)
nspp=ncol(y)
array.lsk.init=res$array.lsk
phi1=matrix(res$phi[nrow(res$phi),],ncomm.init,ncol(y))

#look at convergence
plot(res$llk,type='l')
seq.conv=350:length(res$llk)
plot(res$llk[seq.conv],type='l')

#determine optimal number of groups
nlk=apply(array.lsk.init,c(1,3),sum)
theta1=nlk/apply(nlk,1,sum)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(theta1)
ncomm=4
seq1=1:ncomm
sum(nlk[,seq1]/sum(nlk))

#re-distribute individuals within array.lsk.init that are in eliminated communities
array.lsk=array.lsk.init[,,seq1]
for (i in 1:nloc){
  for (j in 1:nspp){
    tmp=array.lsk.init[i,j,-seq1]
    n=sum(tmp)
    if (n>0){
      prob=theta1[i,seq1]*phi1[seq1,j]
      prob=prob/sum(prob)
      z=rmultinom(1,size=n,prob=prob)
      array.lsk[i,j,]=array.lsk[i,j,]+z
    }
  }
}

#export array.lsk
setwd('U:\\GIT_models\\git_LDA_MS')
dat1=matrix(array.lsk,nloc*nspp*ncomm,1)
write.csv(dat1,'array lsk.csv',row.names=F)

#export phi
tmp=matrix(1:(ncomm.init*nspp),ncomm.init,nspp)
ind1=tmp[-seq1,] #indicators for superfluous groups
write.csv(res$phi[seq.conv,-ind1],'phi step1.csv',row.names=F)

#------------------------------------------------------------------
#are the estimate phi's good?
ordem=numeric()
for (i in 1:ncomm){
  tmp=rep(0,ncomm)
  for (j in 1:ncomm){
    tmp[j]=cor(phi.true[i,],phi1[j,])
  }
  ind=which(tmp==max(tmp))
  ordem=c(ordem,ind)
}
rango=range(c(phi.true,phi1[ordem,]))
plot(phi.true,phi1[ordem,],xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#are the estimate theta's good?
rango=range(c(theta.true,theta1))
plot(theta.true,theta1[,ordem],xlim=rango,ylim=rango)
lines(rango,rango,col='red')