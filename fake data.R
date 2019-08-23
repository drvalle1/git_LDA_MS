rm(list=ls(all=TRUE))
library(MCMCpack)
set.seed(5)

nloc=500
nspp=40
ncommun=5

#design matrix
xmat=matrix(rnorm(nloc*ncommun,sd=2),nloc,ncommun)

#parameters
lambda.true=lambda=runif(ncommun,min=10,max=12)
betas.true=betas=diag(1,ncommun)

#get means
lambda1=matrix(lambda,nloc,ncommun,byrow=T)
media=exp(log(lambda1)+xmat%*%betas); range(media)

#generate N_lk
nlk=matrix(NA,nloc,ncommun)
for (i in 1:ncommun){
  nlk[,i]=rpois(nloc,media[,i])
}
nlk.true=nlk; boxplot(nlk)
z=nlk/apply(nlk,1,sum); apply(z,1,sum); boxplot(z); apply(z,2,range)
nl=apply(nlk,1,sum)
hist(nl)
sum(nl)

#generate phi (assuming that each species is strongly present in a single group) 
phi=rdirichlet(ncommun,alpha=rep(0.1,nspp))

# for (i in 1:nspp){ #add some zeroes
#   ind=sample(1:ncommun,size=1)
#   tmp[ind,i]=runif(1,min=0.5,max=1)
# }
# phi=tmp/matrix(rowSums(tmp),ncommun,nspp) #re-scale to make sure it sums to 1
# round(phi[,1:20],2)
# table(round(phi,2))

unique(rowSums(phi))
phi.true=phi
image(phi)

#generate actual observations y
y=matrix(NA,nloc,nspp)
nks=matrix(0,ncommun,nspp)
for (i in 1:nloc){
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
phi.estim=nks/matrix(rowSums(nks),ncommun,nspp,)
plot(phi.true,phi.estim)
nks.true=nks

#export results
setwd('U:\\GIT_models\\git_LDA_MS')
nome=paste('fake data',ncommun,'.csv',sep='')    
colnames(y)=paste('spp',1:nspp,sep='')
rownames(y)=paste('loc',1:nloc,sep='')
write.csv(y,nome,row.names=F)    

nome=paste('fake data xmat',ncommun,'.csv',sep='')    
write.csv(xmat,nome,row.names=F)    