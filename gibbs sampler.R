# rm(list=ls(all=TRUE))
library('Rcpp')
library('RcppArmadillo')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data5.csv',as.is=T)
xmat=data.matrix(read.csv('fake data xmat5.csv',as.is=T))
y=data.matrix(dat)

#basic settings
ncomm=5
ngibbs=1000
nburn=ngibbs/2
nparam=ncol(xmat)
nloc=nrow(dat)
nspp=ncol(dat)

#priors
phi.prior=0.1
lambda.a=lambda.b=0.1
betas=matrix(0,nparam,ncomm)

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
phi.out=matrix(NA,ngibbs,nspp*ncomm)
nlk.out=matrix(NA,ngibbs,nloc*ncomm)
llk.out=rep(NA,ngibbs)
betas.out=matrix(NA,ngibbs,nparam*ncomm)

#useful stuff for MH algorithm
accept1=list(betas=matrix(0,nparam,ncomm))
jump1=list(betas=matrix(1,nparam,ncomm))
accept.output=50
nadapt=ngibbs/2

#run gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #get log mean
  llambda=matrix(log(lambda),nloc,ncomm,byrow=T)
  lmedia=llambda+xmat%*%betas
  
  #sample z
  tmp=samplez.R(lphi=log(phi), lmedia=lmedia, array.lsk=array.lsk, y=y,
              nlk=nlk, ncomm=ncomm,nloc=nloc, nspp=nspp)
  nlk=tmp$nlk
  array.lsk=tmp$array.lsk
  nks=t(apply(array.lsk,2:3,sum))
    
  #sample phi
  # phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
  phi=phi.true
  
  #sample betas
  tmp=sample.betas(lambda.a=lambda.a,lambda.b=lambda.b,nlk=nlk,xmat=xmat,betas=betas,
                   ncomm=ncomm,nparam=nparam,jump1=jump1$betas)
  betas=tmp$betas
  accept1$betas=accept1$betas+tmp$accept
  # betas=betas.true

  #sample lambdas
  lambda=sample.lambdas(lambda.a=lambda.a,lambda.b=lambda.b,nlk=nlk,ncomm=ncomm,nloc=nloc,
                        xmat=xmat,betas=betas)
  # lambda=lambda.true

  #adaptive MH
  if (i%%accept.output==0 & i<nadapt){
    k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
    accept1=k$accept1
    jump1=k$jump1
  }
  
  #calculate loglikelihood
  p1=dpois(nlk,matrix(lambda,nloc,ncomm,byrow=T),log=T)
  phi.tmp=phi; phi.tmp[phi.tmp<0.00001]=0.00001
  p2=nks*log(phi.tmp)
  
  #store results  
  llk.out[i]=sum(p1)+sum(p2)
  phi.out[i,]=phi
  lambda.out[i,]=lambda
  nlk.out[i,]=nlk
  betas.out[i,]=betas
}
