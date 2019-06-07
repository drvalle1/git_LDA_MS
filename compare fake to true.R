tmp=res$theta[nrow(res$theta),]
theta=matrix(tmp,nloc,ncomm)
boxplot(theta)

plot(NA,xlim=c(0,nloc),ylim=c(0,1))
cores=c('black','red','blue','green','cyan','grey','pink','yellow','brown','orange')
for (i in 1:ncomm){
  lines(1:nloc,theta[,i],col=cores[i])
}

tmp=res$psi[nrow(res$psi),]
psi=matrix(tmp,nloc,ncomm-1)
boxplot(psi)
abline(h=mu0,col='red')

mu=res$mu[nrow(res$mu),]
plot(mu,type='h')