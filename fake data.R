rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

setwd('U:\\GIT_models\\git_LDA_MS')
sourceCpp('aux1.cpp')

nloc=1000
nspp=100
ncommun=4
mu=0
sd1=0.1
sig2=sd1^2  

#generate covariates
tmp1=c(seq(from=2,to=0,length.out=nloc/2),rep(0,nloc*1/2))
tmp2=c(rep(0,nloc/8),seq(from=0,to=2,length.out=nloc/8),seq(from=2,to=0,length.out=nloc/8),rep(0,nloc*5/8))
tmp3=c(rep(0,nloc/4),seq(from=0,to=2,length.out=nloc/8),seq(from=2,to=0,length.out=nloc/8),rep(0,nloc*4/8))
xmat=cbind(tmp1,tmp2,tmp3)
colnames(xmat)=paste('cov',1:(ncommun-1),sep='')

#look at xmat
plot(NA,NA,ylim=range(xmat),xlim=c(1,nloc))
for (i in 1:ncol(xmat)) lines(1:nloc,xmat[,i],col=i)

#generate betas
seq1=c(2,3,2)
betas=diag(seq1,ncommun-1)
betas[betas==0]=-3
betas.true=betas

#generate deltas
tmp=rnorm(nloc*(ncommun-1),mean=mu+xmat%*%betas,sd=sd1)
delta.true=delta=matrix(tmp,nloc,ncommun-1)
range(delta) #this is important: defines truncation points

#generate thetas
prob=1/(1+exp(-delta))
vmat=cbind(prob,1)
theta.true=theta=convertVtoTheta(vmat,rep(1,nloc))
boxplot(theta) #this is important: it determines if the model will be able to disentagle these groups

#create true curves
par(mfrow=c(2,2),mar=rep(2,4))
rango=range(xmat)
nsim=1000
ncov.val=500
seq1=seq(from=rango[1],to=rango[2],length.out=ncov.val)
for (i in 1:(ncommun-1)){
  x=matrix(0,ncov.val,ncommun-1)
  x[,i]=seq1
  media=mu+x%*%betas
  res=matrix(NA,ncov.val,nsim)
  for (j in 1:nsim){
    tmp=rnorm(ncov.val*(ncommun-1),mean=media,sd=sd1)
    delta=matrix(tmp,ncov.val,ncommun-1)
    prob=1/(1+exp(-delta))
    vmat=cbind(prob,1)
    tmp=convertVtoTheta(vmat,rep(1,nloc))
    res[,j]=tmp[,i]
  }
  res1=apply(res,1,quantile,c(0.025,0.5,0.975))
  plot(seq1,res1[2,],ylim=c(0,1))
  lines(seq1,res1[1,],lty=3)
  lines(seq1,res1[3,],lty=3)
}

#generate phi  
tmp=matrix(rnorm(ncommun*nspp,mean=0,sd=2),ncommun,nspp)
tmp[tmp<0.1]=0.1
tmp[,1:(2*ncommun)]=cbind(diag(8,ncommun),diag(8,ncommun))
phi=tmp/matrix(rowSums(tmp),ncommun,nspp)
round(phi[,1:20],2)
table(round(phi,2))

unique(rowSums(phi))
phi.true=phi

#generate actual observations y
nl=floor(runif(nloc,min=100,max=200))
nlk=matrix(NA,nloc,ncommun)
nks=matrix(0,ncommun,nspp)
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  nlk[i,]=rmultinom(1,size=nl[i],prob=theta[i,])
  tmp1=rep(0,ncommun)
  for (k in 1:ncommun){
    tmp=rmultinom(1,size=nlk[i,k],prob=phi[k,])
    nks[k,]=nks[k,]+tmp
    tmp1=tmp1+tmp
  }
  y[i,]=tmp1
}
image(y)

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
nome=paste(c('fake data','fake data cov'),ncommun,'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome[1])
write.csv(xmat,nome[2],row.names=F)    

