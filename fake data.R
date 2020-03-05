rm(list=ls(all=TRUE))
library(MCMCpack)
set.seed(13)

nloc=3000
nspp=150
ncommun=5

#design matrix
xmat=matrix(runif(nloc*ncommun,min=-1,max=3),nloc,ncommun)

#pure sites
tmp=matrix(-3,ncommun,ncommun)
diag(tmp)=2
num1=floor(nloc/ncommun)
for (i in 1:300){
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
NBN=10 #when this is large, we get into areas with relatively flat loglikel, giving trouble to the slice sampler
for (i in 1:ncommun){
  nlk[,i]=rnbinom(nloc,mu=media[,i],size=NBN)
}
nlk.true=nlk; boxplot(nlk)
soma=apply(nlk,1,sum); sum(soma==0)
theta.true=theta=nlk/soma; apply(theta,1,sum); 
boxplot(theta); apply(theta,2,range); apply(theta>0.95,2,mean)
nl=apply(nlk,1,sum)
hist(nl)
sum(nl)

plot(media,nlk)

#generate phi:
#- assume that each species is strongly present in a single group
#- Avoid very rare species (species that are almost never present)
tmp=matrix(0,ncommun,nspp)
base=nspp/ncommun
margin1=floor(base*0.2)
seq1=c(seq(from=1,to=nspp,by=base),nspp-margin1)
for (i in 1:(length(seq1)-1)){ #add some zeroes
  seq2=seq1[i]:(seq1[i+1]+margin1)
  tmp[i,seq2]=1
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
phi.estim=nks/rowSums(nks)
rango=range(c(phi,phi.estim))
plot(phi,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango,col='red')

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