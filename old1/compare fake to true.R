betas.estim=matrix(apply(res$betas,2,mean),ncol(dmat),ncomm-1)
sig2.estim=sig2=apply(res$sig2,2,mean)
psi=matrix(apply(res$psi,2,mean),nloc,ncomm-1)
phi.estim=matrix(apply(res$phi,2,mean),ncomm,ncol(dat))

tmp=res$theta[nrow(res$theta),]
theta.estim=matrix(tmp,nloc,ncomm)
boxplot(theta.estim)

#compare thetas
ind1=numeric()
for (i in 1:ncol(theta.true)){
  cor1=numeric()
  for (j in 1:ncol(theta.estim)){
    tmp=cbind(theta.estim[,j],theta.true[,i])
    cor1=c(cor1,cor(tmp)[1,2])
  }
  ind=which(cor1==max(cor1))
  ind1=c(ind1,ind)
}
ind1
theta.estim1=theta.estim[,ind1]

rango=c(0,1)
plot(theta.true,theta.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red',lwd=2)

#look at phi
phi.estim1=phi.estim[ind1,]
rango=range(c(phi.true,phi.estim1))
plot(phi.true,phi.estim1,xlim=rango,ylim=rango)
lines(rango,rango,col='red',lwd=2)

#------------------------------------------------
resp.curve=function(medias,sig2,nsim){
  nobs=nrow(medias)
  ngroups=ncol(medias)+1
  theta=matrix(NA,nobs,ngroups)
  sig.mat=matrix(sqrt(sig2),nsim,ngroups-1,byrow=T)
  for (i in 1:nobs){
    medias.mat=matrix(medias[i,],nsim,ngroups-1,byrow=T)
    tmp=rnorm(nsim*(ngroups-1),mean=medias.mat,sd=sig.mat)
    tmp1=matrix(tmp,nsim,ngroups-1)
    v=pnorm(tmp1)
    res=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
    theta[i,]=apply(res,2,mean)
  }
  theta
}
#------------------------------------------------
#look at response curves
nparam=ncol(dmat)
seq1=seq(from=0,to=6,length.out=100)
for (i in 1:nparam){
  xmat=matrix(6,100,nparam)
  xmat[,i]=seq1
  
  #calculate implied true theta
  media.true=xmat%*%betas.true
  theta=resp.curve(medias=media.true,sig2=sig2.true,nsim=1000)
  plot(NA,NA,xlim=range(seq1),ylim=c(0,1),main=i)
  for (j in 1:ncol(theta)){
    lines(seq1,theta[,j],col=j)
  }
  
  #calculate implied estimated theta
  media.estim=xmat%*%betas.estim
  theta1=resp.curve(medias=media.estim,sig2=sig2.estim,nsim=1000)[,ind1]
  for (j in 1:ncol(theta1)){
    lines(seq1,theta1[,j],col=j,lty=2)
  }
}