rm(list=ls(all=TRUE))
set.seed(1)

#settings
ncomm=10
nsim=20
nloc=1000
lambda=0.1
tmp=exp(-lambda*1:ncomm)
tmp=tmp-min(tmp) #from 0
tmp=tmp/max(tmp) #to 1
muk=20*tmp
plot(1:ncomm,muk)
sd=0.1

#run simulations
for (i in 1:nsim){
  b0k=rnorm(ncomm-1,mean=muk,sd)
  lomega=matrix(NA,nloc,ncomm)
  lomega[,ncomm]=1
  for (j in 1:(ncomm-1)){
    lomega[,j]=rnorm(nloc,mean=b0k[j],sd=1)
  }
  omega=exp(lomega)
  soma=matrix(apply(omega,1,sum),nloc,ncomm)
  theta=omega/soma
  boxplot(theta,main=i)
}
