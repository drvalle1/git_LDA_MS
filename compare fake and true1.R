#look at convergence
plot(res$llk,type='l',ylim=range(res$llk,na.rm=T))

#look at theta
nloc=nrow(y)
theta=matrix(res$theta[nrow(res$theta),],nloc,ncomm)

par(mfrow=c(1,1))
boxplot(theta)

ind=c(1,2,6,9)
theta1=theta[,ind]
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncol(theta1)){
  lines(theta1[,i],col=i)
}

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
par(mfrow=c(2,2))
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
