rm(list=ls(all=TRUE))

ncomm=4
betas=data.frame(b0=c(1,1,0),b1=c(0.5,0.7,1))
nobs=1000
seq1=seq(from=-1,to=5,length.out=nobs)
xmat=cbind(1,seq1)

lomega=matrix(NA,nobs,ncomm-1)
for (i in 1:(ncomm-1)){
  betas1=as.numeric(betas[i,])
  lomega[,i]=xmat%*%betas1
}
omega=cbind(exp(lomega),1)
soma=rowSums(omega)
theta=omega/matrix(soma,nobs,ncomm)

plot(NA,NA,xlim=c(0,nobs),ylim=c(0,1))
for (i in 1:ncomm){
  lines(1:nobs,theta[,i],col=i)
  text(nobs,theta[nobs,i],col=i,i)
}

#interpretation of coefficients is difficult but response curves will generally be meaningful
#comm 1: should increase with x but actually decreases
#comm 2: should increase with x but actually increases and then decreases

#b0j and b1j have to do with the log-odds of group j relative to the baseline group
ratio=theta[,1]/theta[,4]
eqt=exp(betas$b0[1]+betas$b1[1]*seq1)
plot(1:nobs,ratio,type='l')
lines(1:nobs,eqt,col='red')