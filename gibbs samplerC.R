# rm(list=ls(all=TRUE))
library('Rcpp')
library('mvtnorm')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
sourceCpp('aux1.cpp')
ncomm=10

#get data
dat=read.csv('fake data4.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
nspp=ncol(y)
nloc=nrow(y)

#get covariates
xmat=data.matrix(read.csv('fake data cov4.csv',as.is=T))
npar=ncol(xmat)

#useful stuff
hi=0.999999
lo=0.000001

#initial values of parameters
betas=matrix(0,npar,ncomm-1)
theta=matrix(1/ncomm,nloc,ncomm)
phi=matrix(1/nspp,ncomm,nspp)

#MH stuff
accept1=list(betas=matrix(0,npar,ncomm-1))
jump1=list(betas=matrix(1,npar,ncomm-1))
accept.output=100

#gibbs details
ngibbs=1000
theta.out=matrix(NA,ngibbs,ncomm*nloc)
phi.out=matrix(NA,ngibbs,ncomm*nspp)
betas.out=matrix(NA,ngibbs,npar*(ncomm-1))

llk=rep(NA,ngibbs)
options(warn=2)
for (i in 1:ngibbs){
  print(i)   
  
  #sample betas
  tmp=get.betas(xmat=xmat,y=y,
                betas=betas,phi=phi,
                npar=npar,ncomm=ncomm,
                jump=jump1$betas,lo=lo)
  betas=tmp$betas
  theta=tmp$theta
  accept1$betas=accept1$betas+tmp$accept
  
  #sample z
  tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
  nlk=tmp$nlk
  # nlk=nlk.true
  nks=tmp$nks
  # nks=nks.true
  
  #sample phi
  phi=rdirichlet1(alpha=nks+1,ncomm=ncomm,nspp=nspp) 
  # phi[phi>hi]=hi; phi[phi<lo]=lo
  # phi[4,]=phi.true[4,]
  
  #adapt MH
  if (i%%accept.output==0 & i<1000){
    tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output) 
    jump1=tmp$jump1
    accept1=tmp$accept1
  }
  
  #calculate loglikelihood
  prob=theta%*%phi; prob[prob<lo]=lo #to avoid numerical issues

  #store results  
  llk[i]=sum(y*log(prob))
  theta.out[i,]=theta
  phi.out[i,]=phi
  betas.out[i,]=betas
}

plot(llk,type='l',ylim=range(llk,na.rm=T))
