#look at convergence
plot(res$llk,type='l',ylim=range(res$llk,na.rm=T))

#look at theta
nloc=nrow(y)
theta=matrix(res$theta[nrow(res$theta),],nloc,ncomm)

plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncomm){
  lines(theta[,i],col=i)
}

boxplot(theta)

#create true curves
# ind=c(1,3:5)
# betas1=betas[,ind]
# 
# #create true curves
# par(mfrow=c(2,2))
# rango=range(xmat)
# nsim=1000
# ncov.val=500
# seq1=seq(from=rango[1],to=rango[2],length.out=ncov.val)
# for (i in 1:(ncommun-1)){
#   x=cbind(1,matrix(0,ncov.val,ncommun-1))
#   x[,i+1]=seq1
#   media.estim=x%*%betas1
#   media=x%*%betas.true
#   res.estim=res=matrix(NA,ncov.val,nsim)
#   for (j in 1:nsim){
#     #true
#     tmp=rnorm(ncov.val*(ncommun-1),mean=media,sd=sd1)
#     delta=matrix(tmp,ncov.val,ncommun-1)
#     prob=1/(1+exp(-delta))
#     vmat=cbind(prob,1)
#     tmp=convertVtoTheta(vmat,rep(1,nloc))
#     res[,j]=tmp[,i]
#     
#     #estim
#     tmp=rnorm(ncov.val*(ncommun-1),mean=media.estim,sd=sd1)
#     delta=matrix(tmp,ncov.val,ncommun-1)
#     prob=1/(1+exp(-delta))
#     vmat=cbind(prob,1)
#     tmp=convertVtoTheta(vmat,rep(1,nloc))
#     res.estim[,j]=tmp[,i]
#   }
#   res1=apply(res,1,quantile,c(0.025,0.5,0.975))
#   res1.estim=apply(res.estim,1,quantile,c(0.025,0.5,0.975))
#   
#   #true
#   plot(seq1,res1[2,],ylim=c(0,1))
#   lines(seq1,res1[1,],lty=3)
#   lines(seq1,res1[3,],lty=3)
#   
#   #estim
#   lines(seq1,res1.estim[2,],col='red')
#   lines(seq1,res1.estim[3,],lty=3,col='red')
#   lines(seq1,res1.estim[1,],lty=3,col='red')
# }
