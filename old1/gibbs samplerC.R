# rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
dat=read.csv('fake data5.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
nspp=ncol(y)
nloc=nrow(y)

#priors
mu=5
sd1=3.286
sig2=sd1^2  

#useful stuff
ncomm=10
hi=0.999999
lo=0.000001

#initial values of parameters
theta=matrix(1/ncomm,nloc,ncomm)
delta=matrix(0,nloc,ncomm-1)
phi=matrix(1/nspp,ncomm,nspp)
gamma=0.1; 

#MH stuff
accept1=list(delta=matrix(0,nloc,ncomm-1))
jump1=list(delta=matrix(0.3,nloc,ncomm-1))
accept.output=100

#gibbs details
ngibbs=1000
theta.out=matrix(NA,ngibbs,ncomm*nloc)
phi.out=matrix(NA,ngibbs,ncomm*nspp)
llk=rep(NA,ngibbs)
options(warn=2)
for (i in 1:ngibbs){
  print(i)   
  
  #sample z
  tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
  nlk=tmp$nlk
  # nlk=nlk.true
  nks=tmp$nks
  # nks=nks.true
  
  #get parameters
  tmp=get.delta.theta(nlk=nlk,gamma=gamma,ncomm=ncomm,nloc=nloc,
                      delta=delta,sig2=sig2,mu=mu,jump=jump1$delta)
  delta=tmp$delta
  theta=tmp$theta
  accept1$delta=accept1$delta+tmp$accept
  # theta[theta>hi]=hi; theta[theta<lo]=lo
  # theta=theta.true
  
  phi=rdirichlet1(alpha=nks+1,ncomm=ncomm,nspp=nspp) 
  # phi[phi>hi]=hi; phi[phi<lo]=lo
  # phi=phi.true
  
  #adapt MH
  if (i%%accept.output==0 & i<1000){
    tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output) 
    jump1=tmp$jump1
    accept1=tmp$accept1
  }
  
  #calculate loglikelihood
  prob=theta%*%phi
  prob[prob>hi]=hi; prob[prob<lo]=lo

  #store results  
  llk[i]=sum(y*log(prob))
  theta.out[i,]=theta
  phi.out[i,]=phi
}

plot(llk,type='l',ylim=range(llk,na.rm=T))
