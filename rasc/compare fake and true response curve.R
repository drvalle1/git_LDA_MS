#get betas
nloc=nrow(xmat)
nparam=ncol(xmat)
tmp=res$betas[nrow(res$betas),]
betas=matrix(tmp,nparam,(ncomm-1))

#look at response curves
rango=apply(xmat.centered,2,range)
for (i in 2:nparam){
  xmat1=matrix(0,nloc,nparam); xmat1[,1]=1
  xmat1[,i]=seq(from=rango[1,i],to=rango[2,i],length.out=nloc)
  media.true=xmat1%*%betas.true
  media=xmat1%*%betas
  
  vmat.t=cbind(pnorm(media.true),1)
  theta.t=convertVtoTheta(vmat.t,rep(1,nloc))

  vmat.e=cbind(pnorm(media),1)
  theta.e=convertVtoTheta(vmat.e,rep(1,nloc))
  
  par(mfrow=c(2,1))
  plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1),main=paste('true',i))
  for (j in 1:ncol(theta.t)) lines(1:nloc,theta.t[,j],col=j)
  plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1),main=paste('estimated',i))
  for (j in 1:ncol(theta.e)) lines(1:nloc,theta.e[,j],col=j)
  
  par(mfrow=c(1,1))
  ind=which(theta.t[800,]==max(theta.t[800,]))
  plot(1:nloc,theta.t[,ind],main=i,type='l')
  ind=which(theta.e[800,]==max(theta.e[800,]))
  lines(1:nloc,theta.e[,ind],col='red')
  
}