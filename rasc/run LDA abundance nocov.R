rm(list=ls(all=TRUE))
library('Rcpp')
library('RcppArmadillo')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
source('LDA.abundance main function.R')
sourceCpp('aux1.cpp')

#get data
# dat=read.csv('fake data5.csv',as.is=T)
# ind=which(colnames(dat)=='X')
# y=data.matrix(dat[,-ind]); dim(y)

setwd('U:\\GIT_models\\LdaPoisson_nocov')
dat=read.csv('fake data.csv',as.is=T)
y=data.matrix(dat)

ncomm=30
ngibbs=1000
nburn=ngibbs/2
psi=0.01
gamma=0.1
res=LDA.abundance(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,psi=psi,gamma=gamma)

plot(res$llk,type='l')
seq1=250:nrow(res$theta)
theta=colMeans(res$theta[seq1,])
theta1=matrix(theta,nrow(y),ncomm)
boxplot(theta1)

#compare estimated with true
tmp=apply(res$array.lsk,c(1,3),sum)
ncomm=8
nlk.estim=tmp[,1:ncomm]
ordem=numeric()
for (i in 1:ncomm){
  tmp=rep(NA,ncomm)
  for (j in 1:ncomm){
    tmp[j]=cor(cbind(nlk.true[,i],nlk.estim[,j]))[1,2]
  }
  ind=which(tmp==max(tmp))
  ordem=c(ordem,ind)
}

#compare nlk
rango=range(c(nlk.true,nlk.estim))
plot(nlk.true,nlk.estim[,ordem],xlim=rango,ylim=rango)
lines(rango,rango)

#compare phi
tmp=matrix(res$phi[ngibbs/2,],30,ncol(y))
phi.estim=tmp[1:ncomm,]
rango=range(c(phi.true,phi.estim[ordem,]))
plot(phi.true,phi.estim[ordem,],xlim=rango,ylim=rango)
lines(rango,rango)

