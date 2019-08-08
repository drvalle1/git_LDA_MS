rm(list=ls(all=TRUE))
library('Rcpp')
library('RcppArmadillo')
library(inline)
set.seed(4)

#get functions
setwd('U:\\independent studies\\LDA explorations')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data8.csv',as.is=T)
nloc=nrow(dat)

#basic settings
ncomm=8
ngibbs=1000
nburn=ngibbs/2

#priors
phi.prior=0.1
lambda.a=lambda.b=1

#useful stuff
y=data.matrix(dat)
nloc=nrow(y)
nspp=ncol(y)

#initial values
array.lsk=array(0,dim=c(nloc,nspp,ncomm))
for (i in 1:nloc){
  for (j in 1:nspp){
    if (y[i,j]!=0){
      array.lsk[i,j,]=rmultinom(1,size=y[i,j],prob=rep(1/ncomm,ncomm))  
    }
  }
}

#basic test
# z=apply(array.lsk,1:2,sum)
# unique(y-z)

nlk=apply(array.lsk,c(1,3),sum)
nks=t(apply(array.lsk,2:3,sum))

lambda=apply(nlk,2,mean)
phi=matrix(1/nspp,ncomm,nspp)

#to store outcomes from gibbs sampler
lambda.out=matrix(NA,ngibbs,ncomm)
phi.out=matrix(NA,ngibbs,nloc*nspp)
llk=rep(NA,ngibbs)
  
#run gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #get log mean
  llambda=matrix(log(lambda),nloc,ncomm,byrow=T)
  lmedia=llambda
  
  #sample z
  tmp=samplez.R(lphi=log(phi), lmedia=lmedia, array.lsk=array.lsk, y=y,
              nlk=nlk, ncomm=ncomm,nloc=nloc, nspp=nspp)
  nlk=tmp$nlk
  array.lsk=tmp$array.lsk
  nks=t(apply(array.lsk,2:3,sum))
    
  #sample phi
  phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
  
  #sample lambdas
  lambda=sample.lambdas(lambda.a=lambda.a,lambda.b=lambda.b,nlk=nlk,ncomm=ncomm,nloc=nloc)
  
  #calculate loglikelihood
  p1=dpois(nlk,matrix(lambda,nloc,ncomm,byrow=T),log=T)
  p2=nks*log(phi)
  
  #store results  
  llk[i]=sum(p1)+sum(p2)
  phi.out[i,]=phi
  lambda.out[i,]=lambda
}
