rm(list=ls(all=TRUE))
library(MCMCpack)
set.seed(5)

nloc=5000
nspp=100
ncommun=4
base=floor(nloc/(ncommun-2))
nparam=ncommun-1

#we have to make sure that each community dominates at least some locations
xmat=cbind(c(rep(0,nloc*0.2),rep(6,nloc*0.6),runif(nloc*0.2,min=0,max=6)),
           c(rep(6,nloc*0.2),rep(0,nloc*0.2),rep(6,nloc*0.4),runif(nloc*0.2,min=0,max=6)),
           c(rep(6,nloc*0.4),rep(0,nloc*0.2),rep(6,nloc*0.2),runif(nloc*0.2,min=0,max=6)))
betas=matrix(0,nparam,ncommun-1)
diag(betas)=c(-1,-1,-1)
betas.true=betas
mu=3
medias.true=medias=mu+xmat%*%betas
psi=matrix(NA,nloc,ncommun-1)
sig2.true=sig2=c(1,2,1)
for (i in 1:(ncommun-1)){
  psi[,i]=rnorm(nloc,mean=medias[,i],sd=sqrt(sig2[i]))
}
psi.true=psi
apply(psi,2,range)

#generate thetas
theta=matrix(NA,nloc,ncommun)
aux=rep(1,nloc)
for (i in 1:(ncommun-1)){
  tmp=pnorm(psi[,i])
  theta[,i]=tmp*aux
  aux=aux*(1-tmp)
}
theta[,ncommun]=aux
theta.true=theta
  
boxplot(theta)
apply(theta>0.9,2,mean)

unique(apply(theta,1,sum))
for (i in 1:(ncommun-1)) plot(xmat[,i],theta[,i],main=i)

#generate phi  
tmp=rdirichlet(ncommun,alpha=rep(1,nspp))
for (i in 1:nspp){
  ind=sample(1:ncommun,size=ncommun/2)
  tmp[ind,i]=0
}
phi=tmp/matrix(rowSums(tmp),ncommun,nspp)
round(phi[,1:20],2)
table(round(phi,2))

unique(rowSums(phi))
phi.true=phi
image(phi)

#generate actual observations y
nl=floor(runif(nloc,min=200,max=400))
nlk=matrix(NA,nloc,ncommun)
nks=matrix(0,ncommun,nspp)
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  nlk[i,]=rmultinom(1,size=nl[i],prob=theta[i,])
  tmp1=rep(0,nspp)
  for (k in 1:ncommun){
    tmp=rmultinom(1,size=nlk[i,k],prob=phi[k,])
    nks[k,]=nks[k,]+tmp
    tmp1=tmp1+tmp
  }
  y[i,]=tmp1
}
# image(y)

#look at stuff to make sure it makes sense
theta.estim=nlk/matrix(nl,nloc,ncommun)
plot(NA,NA,xlim=c(1,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta.estim[,i],col=i)
nlk.true=nlk

phi.estim=nks/matrix(rowSums(nks),ncommun,nspp,)
plot(phi.true,phi.estim)
nks.true=nks

#export results
setwd('U:\\GIT_models\\git_LDA_MS')
nome=paste('fake data',ncommun,c('',' dmat'),'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome[1])    
write.csv(xmat,nome[2],row.names=F)