#look at convergence
plot(res$llk,type='l',ylim=range(res$llk,na.rm=T))

#look at theta
nloc=nrow(y)
theta=matrix(res$theta[nrow(res$theta),],nloc,ncomm)
boxplot(theta)

#re-order theta
true.ncomm=9
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

#look at sig2
sig2=res$sig2
for (i in 1:ncol(sig2)) {
  plot(sig2[,i],type='l',main=i)
  abline(h=0.1,col='red')
}
sig2=res$sig2[nrow(res$sig2),]
round(sig2,3)
plot(1:length(sig2),sig2)
abline(h=0.1,col='red')

#get parameters
npar=ncol(xmat1)
betas=matrix(res$betas[nrow(res$betas),],npar,ncomm-1)

#create design matrix of interest
seq1=seq(from=-1,to=1,length.out=nloc)
knots=seq(from=-0.5,to=0.5,by=0.5)
mat=seq1
for (i in 1:length(knots)){
  tmp=ifelse(seq1-knots[i]<0,0,seq1-knots[i])
  mat=cbind(mat,tmp)
}
mat.focus=mat

#create design matrix to fill in
seq1=rep(0,nloc)
mat=seq1
for (i in 1:length(knots)){
  tmp=ifelse(seq1-knots[i]<0,0,seq1-knots[i])
  mat=cbind(mat,tmp)
}
mat.fill=mat

#create true curves
par(mfrow=c(2,2),mar=c(3,3,1,1))
cov=4

#get design matrix
xmat2=xmat1
xmat2[,-1]=cbind(mat.fill,mat.fill,mat.fill,mat.fill)
nome=paste0('cov',cov)
ind=grep(nome,colnames(xmat2))
xmat2[,ind]=mat.focus
  
#calculate theta
media=xmat2%*%betas
prob=1/(1+exp(-media))
vmat=cbind(prob,1)
theta=convertVtoTheta(vmat,rep(1,nloc))

ind=c(1,2,3,6)  
#plot results
for (i in 1:length(ind)) plot(seq(from=-1,to=1,length.out=nloc),theta[,ind[i]],type='l',main=cov)
