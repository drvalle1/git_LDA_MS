ind=1:3
# ind=c(3,1,2)
betas1=betas[,ind]

par(mfrow=c(2,2))
for (i in 1:3) plot(betas.true[,i],betas1[,i])

table(jump1$delta)
boxplot(jump1$delta)

boxplot(theta)

plot(delta.true,delta[,ind])

#create true curves
par(mfrow=c(2,2))
rango=range(xmat)
nsim=1000
ncov.val=500
seq1=seq(from=rango[1],to=rango[2],length.out=ncov.val)
for (i in 1:(ncommun-1)){
  x=matrix(0,ncov.val,ncommun-1)
  x[,i]=seq1
  media.estim=mu+x%*%betas1
  media=mu+x%*%betas.true
  res.estim=res=matrix(NA,ncov.val,nsim)
  for (j in 1:nsim){
    #true
    tmp=rnorm(ncov.val*(ncommun-1),mean=media,sd=sd1)
    delta=matrix(tmp,ncov.val,ncommun-1)
    prob=1/(1+exp(-delta))
    vmat=cbind(prob,1)
    tmp=convertVtoTheta(vmat,rep(1,nloc))
    res[,j]=tmp[,i]
    
    #estim
    tmp=rnorm(ncov.val*(ncommun-1),mean=media.estim,sd=sd1)
    delta=matrix(tmp,ncov.val,ncommun-1)
    prob=1/(1+exp(-delta))
    vmat=cbind(prob,1)
    tmp=convertVtoTheta(vmat,rep(1,nloc))
    res.estim[,j]=tmp[,i]
  }
  res1=apply(res,1,quantile,c(0.025,0.5,0.975))
  res1.estim=apply(res.estim,1,quantile,c(0.025,0.5,0.975))
  
  #true
  plot(seq1,res1[2,],ylim=range(res1))
  lines(seq1,res1[1,],lty=3)
  lines(seq1,res1[3,],lty=3)
  
  #estim
  lines(seq1,res1.estim[2,],col='red')
  lines(seq1,res1.estim[3,],lty=3,col='red')
  lines(seq1,res1.estim[1,],lty=3,col='red')
}
