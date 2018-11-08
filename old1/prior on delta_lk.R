rm(list=ls(all=TRUE))
gamma=0.1

expected.fun=function(x){
  tmp=1/(1+exp(-x))
  tmp*dnorm(x,mean=mu,sd=sd1)
}
var.fun=function(x){
  tmp=1/(1+exp(-x))
  ((tmp-expected)^2)*dnorm(x,mean=mu,sd=sd1)
}

combo=expand.grid(mu=seq(from=0.01,to=5,length.out=100),
                  sd1=seq(from=0.01,to=5,length.out=100))
combo$expected=combo$var1=NA
for (i in 1:nrow(combo)){
  print(i)
  mu=combo$mu[i]
  sd1=combo$sd1[i]
  z=integrate(expected.fun,lower=-100,upper=100,stop.on.error=F)
  expected=z$value
  tmp=integrate(var.fun,lower=-100,upper=100,stop.on.error=F)
  var1=tmp$value
  if (tmp$message!="the integral is probably divergent") combo$var1[i]=var1
  if (z$message!="the integral is probably divergent")   combo$expected[i]=expected
}
expected=1/(1+gamma)
denom=((2+gamma)^2)*(3+gamma)
var1=(1+gamma)/denom

#eliminate NA's
cond=!is.na(combo$expected) & !is.na(combo$var1)
combo1=combo[cond,]

#calculate differences
combo1$diff.exp=abs(combo1$expected-expected)
combo1$diff.var=abs(combo1$var1-var1)
combo1$soma=combo1$diff.exp+combo1$diff.var

#choose set with close mu's 
ind=which(combo1$soma==min(combo1$soma))
combo1[ind,]

#compare distributions
mu=combo1$mu[ind] #5
sd1=combo1$sd1[ind] #3.286
n=100000
z=rnorm(n,mean=mu,sd=sd1)
z1=1/(1+exp(-z))
mean(z1); var(z1)
expected; var1

plot(density(z1,from=0,to=1),type='l')
x=seq(from=0.0001,to=0.9999,by=0.001)
lines(x,dbeta(x,1,gamma),col='red',lwd=3)