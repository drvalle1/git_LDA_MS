rm(list=ls()) 
nobs=10000
n=100
mu=20
p=n/(n+mu)
x=rnbinom(nobs,size=n,prob=p)

mean(x); mu
var(x); n*(1-p)/(p^2)