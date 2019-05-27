#look at convergence
plot(res$llk,type='l',ylim=range(res$llk,na.rm=T))

#look at theta
nloc=nrow(y)
theta=matrix(res$theta[nrow(res$theta),],nloc,ncomm)
boxplot(theta)

#re-order theta
true.ncomm=5
res1=matrix(NA,true.ncomm,ncomm)
for (i in 1:true.ncomm){
  for (j in 1:ncomm){
    res1[i,j]=cor(cbind(theta.true[,i],theta[,j]))[2,1]
  }
}
ind=numeric()
for (i in 1:true.ncomm) ind=c(ind,which(res1[i,]==max(res1[i,])))
theta1=theta[,ind]

plot(theta.true,theta1)

rango=c(0,1)
plot(NA,NA,xlim=c(0,nrow(theta.true)),ylim=rango)
seq1=1:nrow(theta.true)
for (i in 1:true.ncomm) {
  lines(seq1,theta.true[,i],col=i)
  lines(seq1,theta1[,i],col=i,lty=2,lwd=2) 
}

rango=c(0,1)
for (i in 1:true.ncomm){
  plot(theta.true[,i],theta1[,i],xlim=rango,ylim=rango,main=i)
  lines(rango,rango,col='red',lwd=2)
}
