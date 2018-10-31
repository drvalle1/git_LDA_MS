rm(list=ls(all=TRUE))
gamma=0.1
expected=1/(1+gamma)
denom=((2+gamma)^2)*(3+gamma)
var1=(1+gamma)/denom

dkl.fun=function(x){
  dbeta(x,1,gamma)*(dbeta(x,1,gamma,log=T)-dnorm(x,mean=mu,sd=sd1,log=T))
}

combo=expand.grid(mu=seq(from=0.01,to=5,length.out=100),
                  sd1=seq(from=0.01,to=5,length.out=100))
combo$dkl=NA
for (i in 1:nrow(combo)){
  print(i)
  mu=combo$mu[i]
  sd1=combo$sd1[i]
  z=integrate(dkl.fun,lower=0.0001,upper=0.9999,stop.on.error=F)
  if (z$message!="the integral is probably divergent") combo$dkl[i]=z$value
}

#eliminate NA's
cond=!is.na(combo$dkl)
combo1=combo[cond,]

#calculate differences
ind=which(combo1$dkl==min(combo1$dkl))
combo2=combo1[ind,]
combo2

#compare distributions
mu=combo2$mu #0.867
sd1=combo2$sd1 #0.262
n=100000
z=rnorm(n,mean=mu,sd=sd1)
z1=pnorm(z)
mean(z1); var(z1)
expected; var1

plot(density(z1,from=0,to=1),type='l')
x=seq(from=0.0001,to=0.9999,by=0.001)
lines(x,dbeta(x,1,gamma),col='red',lwd=3)