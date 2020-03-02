rm(list=ls(all=TRUE)) 
set.seed(1)

n=10 #bigger n = smaller variance
mean1=2
p=n/(mean1+n)

n1=10000
y=rnbinom(n1,mu=mean1,size=n)
mean(y)
var(y); n*(1-p)/(p^2)

#look at density
x=0:100
prob1.true=dnbinom(x,mu=mean1,size=n)
tmp=lgamma(x+n)-lgamma(n)-lgamma(x+1)+n*log(p)+x*log(1-p)
prob1.estim=exp(tmp)
plot(x,prob1.true)
points(x,prob1.estim,col='red',cex=0.5)
