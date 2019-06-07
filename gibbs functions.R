get.w=function(nlk,psi,ncomm,nloc){
  soma=matrix(NA,nloc,ncomm-1)
  ge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  pnorm1=pnorm(0,mean=psi,sd=1)
  
  soma=getw(ge=ge, psi=psi, pnorm1=pnorm1,nlk=nlk,nloc=nloc,ncomm=ncomm)
  list(soma=soma,nge=ge)
}    
#-------------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#---------------------------------------------
get.psi=function(soma,nlk,sig2,mu,ncomm,theta,nks,change1,nburn,i){
  #re-order groups from largest to smallest
  if (i%%change1 == 0 & i<nburn){
    med=colMeans(theta)
    ind=order(med,decreasing=T)
    nlk=nlk[,ind]
    nks=nks[ind,]
    ind1=ind[ind!=ncomm]
    soma=soma[,ind1]
    mu=mu[ind1]
  }
  
  nge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  prec=nge[,-ncomm]+(1/sig2)
  var1=1/prec
  mu1=matrix(mu,nloc,ncomm-1,byrow=T)
  pmedia=soma+(mu1/sig2)
  tmp=rnorm(nloc*(ncomm-1),mean=var1*pmedia,sd=sqrt(var1))
  psi=matrix(tmp,nloc,ncomm-1)
  
  list(psi=psi,nlk=nlk,nks=nks,soma=soma,mu=mu)
}
#---------------------------------------------
get.mu=function(psi,nloc,sd0,mu0,sig2,ncomm){
  prec=(nloc/sig2)+(1/(sd0^2))
  var1=1/prec
  soma=colSums(psi)
  pmedia=(soma/sig2)+(mu0/(sd0^2))
  rnorm(ncomm-1,mean=var1*pmedia,sd=sqrt(var1))
}