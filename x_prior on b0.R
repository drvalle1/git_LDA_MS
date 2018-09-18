#LDA MS prior
gamma=0.1
ncomm=10
seq1=1:(ncomm-1)
mu.c=(ncomm-seq1-1)*log(1+gamma)-(ncomm-seq1)*log(gamma)
plot(seq1,mu.c)

tmp=exp(mu.c)
prob.lda.ms=tmp/(1+sum(tmp))
plot(1:9,prob.lda.ms)

#LDA TSB prior
seq1=0:(ncomm-1)
v=1/(1+gamma)
prob.lda.tsb=v*(1-v)^seq1
points(1:ncomm,prob.lda.tsb,col='red',cex=0.6)

#stochastic test
nsim=10000
res=matrix(NA,nsim,ncomm)
for (i in 1:nsim){
  betas=rnorm(ncomm-1,mean=mu.c,sd=1)
  omega=c(exp(betas),1)
  res[i,]=omega/sum(omega)
}
boxplot(res)
seq1=1:(ncomm-1)
points(seq1,prob.lda.ms,col='red',pch=19)