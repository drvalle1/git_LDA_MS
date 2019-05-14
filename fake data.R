rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

setwd('U:\\GIT_models\\git_LDA_MS')
sourceCpp('LDA_MS_c.cpp')

nloc=1000
nspp=200
ncommun=9

#generate covariates
dist=nloc/(ncommun+1)
init=1

xmat=matrix(rnorm(nloc*(ncommun-1)),nloc,ncommun-1)
xmat=cbind(2,xmat)
colnames(xmat)=paste('cov',0:(ncommun-1),sep='')

#look at xmat
plot(NA,NA,ylim=range(xmat),xlim=c(1,nloc),main='covariates')
for (i in 2:ncol(xmat)) lines(1:nloc,xmat[,i],col=i)

#standardize xmat
media1=apply(xmat,2,mean)
sd1=apply(xmat,2,sd)
xmat[,1]=1
for (i in 2:ncol(xmat)){
  xmat[,i]=(xmat[,i]-media1[i])/sd1[i]
}
apply(xmat,2,mean)
apply(xmat,2,sd)
xmat.centered=xmat

#generate betas
betas=matrix(NA,ncol(xmat),ncommun-1)
betas[1,]=-1
betas[-1,]=diag(1,ncol(xmat)-1)
betas=betas*3
betas.true=betas

#generate thetas
vmat.true=vmat=cbind(pnorm(xmat.centered%*%betas.true),1)
theta.true=theta=convertVtoTheta(vmat,rep(1,nloc))
#this is important: it determines if the model will be able to disentagle these groups
boxplot(theta,ylim=c(0,1)) 

#create true curves
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncommun){
  lines(theta[,i],col=i)
}
apply(theta,2,max)

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
nl=floor(runif(nloc,min=100,max=400))
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
nome=paste(c('fake data','fake data cov'),'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome[1])
write.csv(xmat,nome[2],row.names=F)    

