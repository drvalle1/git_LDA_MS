rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
setwd('U:\\independent studies\\LDA_MS\\fake data')
dat=read.csv('fake data 5.csv',as.is=T)

#useful stuff
ind=which(colnames(dat)=='loc.id')
dat1=data.matrix(dat[,-ind])
minus.dat1=1-dat1
# image(data.matrix(dat1))
loc.id=dat$loc.id
nloc=max(loc.id)
xmat=matrix(1,nloc,1)
t.xmat=t(xmat)
xtx=t.xmat%*%xmat
nparam=ncol(xmat)

nspp=ncol(dat1)
ncomm=10
nlinhas=nrow(dat1)
hi=0.999
lo=1-hi
o.lo=log(lo/(1-lo))
o.hi=log(hi/(1-hi))

#get true phi
# setwd('U:\\independent studies\\LDA_MS\\fake data')
# phi.true=read.csv('phi 5.csv',as.is=T)
# phi.true1=rbind(as.matrix(phi.true),0.01,0.01,0.01,0.01,0.01)
# phi.true1[phi.true1>hi]=hi
# phi.true1[phi.true1<lo]=lo

#initial values for parameters
omega=matrix(0,nloc,ncomm)
omega[,ncomm]=o.hi
vmat=exp(omega)/(1+exp(omega))
prod=rep(1,nloc)
theta=convertVtoTheta(vmat,prod)
betas=matrix(0,nparam,ncomm-1)
phi=matrix(0.5,ncomm,nspp)
z=matrix(sample(1:ncomm,size=nlinhas*nspp,replace=T),nlinhas,nspp)
sig2=1

#priors
a.phi=1
b.phi=1
sig2.b=sig2.a=0.1
sig2.a1=(nloc+2*sig2.a)/2

#gibbs details
ngibbs=1000
theta.out=matrix(NA,ngibbs,ncomm*nloc)
phi.out=matrix(NA,ngibbs,ncomm*nspp)
betas.out=matrix(NA,ngibbs,(ncomm-1)*nparam)
llk=rep(NA,ngibbs)

#for MH
jump=list(omega=matrix(0.5,nloc,ncomm))
accept=list(omega=matrix(0,nloc,ncomm))
accept.output=50

options(warn=2)
for (i in 1:ngibbs){
  print(i)   

  #calculate summaries to sample phi
  tmp=getks(z=z, ncommun=ncomm, dat=dat1)
  nks1=tmp$nks1
  nks0=tmp$nks0
  phi=matrix(rbeta(nspp*ncomm,nks1+a.phi,nks0+b.phi),ncomm,nspp) #phi.true#
  # phi=phi.true1
  
  #sample z
  rand.u=matrix(runif(nlinhas*nspp),nlinhas,nspp)
  z=samplez(ltheta=log(theta), l1minustheta=log(1-phi), lphi=log(phi), 
            dat1=dat1, minus_dat1=minus.dat1,
            locid=loc.id,randu=rand.u, ncommun=ncomm, nloc=nloc)

  #sample omega
  nlk=getlk(z=z,locid=loc.id, ncommun=ncomm, nloc=nloc)
  tmp=sample.omega(ncomm=ncomm,nloc=nloc,xmat=xmat,prod=prod,
                   jump=jump$omega,o.lo=o.lo,o.hi=o.hi,
                   omega=omega,betas=betas,nlk=nlk,sig2=sig2)
  omega=tmp$omega
  accept$omega=accept$omega+tmp$accept
  
  #get corresponding theta
  e.omega=exp(omega)
  vmat=e.omega/(1+e.omega)
  theta=convertVtoTheta(vmat,prod)
  
  #sample beta
  betas=rep(3.653,ncomm-1)#sample.betas(ncomm,omega,sig2,xtx,nparam)
    
  #sample sig2
  sig2=2.959^2
  # sig2=sample.sig2(xmat=xmat,betas=betas,omega=omega,sig2.b=sig2.b,a1=sig2.a1)
  
  #sample 
  prob=theta%*%phi
  prob[prob>hi]=hi; prob[prob<lo]=lo
  prob1=prob[loc.id,]

  #adapt MCMC
  if (i<1000 & i%%accept.output==0){
    tmp=print.adapt(accept1=accept,jump1=jump,accept.output=accept.output)
    accept=tmp$accept1
    jump=tmp$jump1
  }
  
  #store results  
  llk[i]=sum(dat1*log(prob1)+minus.dat1*log(1-prob1))
  theta.out[i,]=theta
  phi.out[i,]=phi
  betas.out[i,]=betas
}

plot(llk,type='l',ylim=range(llk,na.rm=T))

