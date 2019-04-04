rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

nloc=1000
nspp=100
ncommun=4

#generate covariates
xmat=cbind(1,
           runif(nloc,min=-1,max=1),
           runif(nloc,min=-1,max=1),
           runif(nloc,min=-1,max=1),
           runif(nloc,min=-1,max=1))

#generate betas
betas=matrix(NA,ncol(xmat),ncommun)
betas[,1]=c(-1,1,0,0,0)*3
betas[,2]=c(-1,0,1,0,0)*3
betas[,3]=c(-1,0,0,1,0)*3
betas[,4]=c(-1,0,0,0,1)*3
betas.true=betas

#generate thetas
tmp=exp(xmat%*%betas)
soma=matrix(apply(tmp,1,sum),nloc,ncommun)
theta=tmp/soma
par(mfrow=c(1,1))
boxplot(theta)

#look at raw response curves
par(mfrow=c(2,2))
for (i in 1:ncommun){
  plot(xmat[,i+1],theta[,i],main=i)  
}

#look at true response curves
par(mfrow=c(2,2),mar=c(3,3,1,1))
for (i in 1:ncommun){
  xmat1=xmat
  xmat1[,2:5]=0
  xmat1[,i+1]=seq(from=-1,to=1,length.out=nloc)
  tmp=exp(xmat1%*%betas)
  soma=matrix(apply(tmp,1,sum),nloc,ncommun)
  theta.hyp=tmp/soma
  for (j in 1:4) plot(xmat1[,i+1],theta.hyp[,j],type='l',main=i)
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
par(mfrow=c(1,1))
image(y)

#look at stuff to make sure it makes sense
nlk.true=nlk

phi.estim=nks/matrix(rowSums(nks),ncommun,nspp)
plot(phi.true,phi.estim)
nks.true=nks

#export results
setwd('U:\\GIT_models\\git_LDA_MS')
nome=paste(c('fake data response curve','fake data response curve cov'),ncommun,'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome[1])
write.csv(xmat,nome[2],row.names=F)    

