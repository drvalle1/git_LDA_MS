rm(list=ls(all=TRUE))
set.seed(4)

nloc=1000
nspp=100
ncommun=5
base=floor(nloc/(ncommun-2))

#generate thetas
x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)
  
init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)
  
theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta
  
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)
  
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
nome=paste('fake data',ncommun,'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome)    
