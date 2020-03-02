rm(list=ls(all=TRUE))
library(MCMCpack)
set.seed(28)

nloc=3000
nspp=100
ncommun=4

#design matrix
xmat=matrix(runif(nloc*ncommun,min=-1,max=3),nloc,ncommun)

#pure sites
tmp=matrix(-3,ncommun,ncommun)
diag(tmp)=1
num1=floor(nloc/ncommun)
for (i in 1:200){
  seq1=(ncommun*(i-1)+1):(ncommun*i)
  xmat[seq1,]=tmp
}
image(xmat)
xmat=cbind(1,xmat)

#parameters
b0=log(runif(ncommun,min=5,max=8))
betas.true=betas=rbind(b0,diag(1,ncommun))

#get means
media.true=media=exp(xmat%*%betas); range(media)
head(media)

#generate N_lk
par(mfrow=c(1,1),mar=rep(4,4))
nlk=matrix(NA,nloc,ncommun)
for (i in 1:ncommun){
  nlk[,i]=rpois(nloc,media[,i])
}
nlk.true=nlk; boxplot(nlk)
soma=apply(nlk,1,sum); sum(soma==0)
z=nlk/soma; apply(z,1,sum); 
boxplot(z); apply(z,2,range); apply(z>0.9,2,mean)
nl=apply(nlk,1,sum)
hist(nl)
sum(nl)

plot(media,nlk)

#generate phi (assuming that each species is strongly present in a single group) 
tmp=rdirichlet(ncommun,rep(0.01,nspp))
apply(tmp,1,sum)

for (i in 1:nspp){ #add some zeroes
  ind=sample(1:ncommun,size=1)
  tmp[ind,i]=runif(1,min=0.1,max=1)
}
phi=tmp/matrix(rowSums(tmp),ncommun,nspp) #re-scale to make sure it sums to 1
phi.true=phi
# round(phi[,1:20],2)
# table(round(phi,2))

unique(rowSums(phi))
phi.true=phi
image(phi)

#generate actual observations y
array.lsk=array(0,dim=c(nloc,nspp,ncommun))
for (i in 1:nloc){
  for (k in 1:ncommun){
    array.lsk[i,,k]=rmultinom(1,size=nlk[i,k],prob=phi[k,])
  }
}
array.lsk.true=array.lsk
y=apply(array.lsk,c(1,2),sum)
nks=t(apply(array.lsk,c(2,3),sum))
image(y)
plot(phi,nks)

#checking if it makes sense
plot(apply(array.lsk,c(1,3),sum),nlk)
lines(c(0,1000),c(0,1000))

#look at stuff to make sure it makes sense
phi.estim=nks/matrix(rowSums(nks),ncommun,nspp)
plot(phi.true,phi.estim)
nks.true=nks

#export results
setwd('U:\\GIT_models\\git_LDA_MS')
nome='fake data.csv'
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome,row.names=F)    

nome='fake data xmat.csv'
write.csv(xmat,nome,row.names=F)    