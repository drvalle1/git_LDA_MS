betas.estim=matrix(apply(res$betas,2,mean),ncol(dmat),ncomm-1)
sig2=apply(res$sig2,2,mean)
psi=matrix(apply(res$psi,2,mean),nloc,ncomm-1)

tmp=res$theta[nrow(res$theta),]
theta=matrix(tmp,nloc,ncomm)
boxplot(theta)

#look at response curves
nparam=ncol(dmat)
seq1=seq(from=0,to=6,length.out=100)
for (i in 1:nparam){
  xmat=matrix(0,100,nparam)
  xmat[,i]=seq1
  
  #calculate implied true theta
  media.true=xmat%*%betas.true
  v=pnorm(media.true) 
  theta=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
  plot(NA,NA,xlim=range(seq1),ylim=c(0,1),main=i)
  for (j in 1:ncomm){
    lines(seq1,theta[,j],col=j)
  }
  
  #calculate implied estimated theta
  media.estim=xmat%*%betas.estim
  v=pnorm(media.estim) 
  theta=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
  for (j in 1:ncomm){
    lines(seq1,theta[,j],col=j,lty=2)
  }
}