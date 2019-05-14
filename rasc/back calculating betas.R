rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

setwd('U:\\GIT_models\\git_LDA_MS')
sourceCpp('LDA_MS_c.cpp')
nloc=1000
nspp=100
ncommun=6

#generate covariates
tmp1=c(seq(from=2,to=0,length.out=nloc/5),rep(0,4*nloc/5))
tmp2=c(seq(from=0,to=2,length.out=nloc/5),seq(from=2,to=0,length.out=nloc/5),rep(0,3*nloc/5))
tmp3=c(rep(0,nloc/5),seq(from=0,to=2,length.out=nloc/5),seq(from=2,to=0,length.out=nloc/5),rep(0,2*nloc/5))
tmp4=c(rep(0,2*nloc/5),seq(from=0,to=2,length.out=nloc/5),seq(from=2,to=0,length.out=nloc/5),rep(0,nloc/5))
tmp5=c(rep(0,3*nloc/5),seq(from=0,to=2,length.out=nloc/5),seq(from=2,to=0,length.out=nloc/5))
xmat=cbind(1,tmp1,tmp2,tmp3,tmp4,tmp5)*2
colnames(xmat)=paste('cov',0:(ncommun-1),sep='')

#look at xmat
plot(NA,NA,ylim=range(xmat),xlim=c(1,nloc),main='covariates')
for (i in 2:ncol(xmat)) lines(1:nloc,xmat[,i],col=i)

#generate betas
betas=matrix(NA,ncol(xmat),ncommun-1)
betas[,1]=c(-1,1,0,0,0,0)*3
betas[,2]=c(-1,0,1,0,0,0)*3
betas[,3]=c(-1,0,0,1,0,0)*3
betas[,4]=c(-1,0,0,0,1,0)*3
betas[,5]=c(-1,0,0,0,0,1)*3
betas.true=betas

#generate probs
probs=media=matrix(NA,nloc,ncommun-1)
for (i in 1:(ncommun-1)){
  media[,i]=xmat%*%betas[,i]
  probs[,i]=1/(1+exp(-media[,i]))
}
probs.true=probs
range(probs)

#can I retrieve betas in a straightforward way if I know probs? solve(xtx)%*%t(x)%*%logit
p1=solve(t(xmat)%*%xmat)
betas.estim=p1%*%t(xmat)%*%media
plot(betas.true,betas.estim)
