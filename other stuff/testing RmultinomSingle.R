rm(list=ls(all=TRUE))
library('Rcpp')
library('RcppArmadillo')
library(inline)
set.seed(4)

#get functions
setwd('U:\\independent studies\\LDA explorations')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

n=1000
runif1=runif(n)
prob=c(0.3,0.1,0.1,0.5)
ncommun=length(prob)
res=rep(0,ncommun)
for (i in 1:n){
  ind=1+RmultinomSingle(runif1=runif1[i], prob=prob, ncommun=ncommun)
  res[ind]=res[ind]+1
}
res/sum(res)

