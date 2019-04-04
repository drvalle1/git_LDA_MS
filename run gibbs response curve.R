# rm(list=ls(all=TRUE))
library('Rcpp')
library('mvtnorm')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_MS')
source('LDA MS functions.R')
source('LDA MS gibbs sampler.R')
sourceCpp('LDA_MS_c.cpp')

#get data
dat=read.csv('fake data response curve4.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)

#get covariates
xmat=data.matrix(read.csv('fake data response curve cov4.csv',as.is=T))
nomes=paste0('cov',1:4)
colnames(xmat)=c('interc',nomes)

#create linear splines for covariates 2:5
xmat1=matrix(1,nrow(xmat),1)
knots=c(-0.5,0,0.5)
for (i in 1:length(nomes)){
  orig=xmat[,nomes[i]]
  res=orig
  for (j in 1:length(knots)){
    tmp=ifelse(orig-knots[j]<0,0,orig-knots[j])
    res=cbind(res,tmp)
  }
  colnames(res)=c(nomes[i],paste0(nomes[i],c('a','b','c')))
  xmat1=cbind(xmat1,res)
}

#run gibbs sampler
ncomm=10
ngibbs=1000
nburnin=ngibbs/2
res=LDA.MS.gibbs(y=y,xmat=xmat1,ncomm=ncomm,ngibbs=ngibbs,nburnin=nburnin)