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
get.psi=function(soma,nlk,sig2,mu,ncomm,theta,nks,change1,nburn,i,dmat,betas){
  #re-order groups from largest to smallest
  # if (i%%change1 == 0 & i<nburn){
  #   med=colMeans(theta)
  #   ind=order(med,decreasing=T)
  #   ind1=ind[ind!=ncomm]
  #   nlk=cbind(nlk[,ind1],nlk[,ncomm])
  #   nks=rbind(nks[ind1,],nks[ncomm,])
  #   soma=soma[,ind1] #I would probably have to re-run getw() to get correct soma. Before this, change psi
  #   sig2=sig2[ind1]
  #   betas=betas[,ind1]
  # }

  medias=mu+dmat%*%betas
  nge=ngreater(nlk=nlk,nloc=nloc,ncommun=ncomm) #greater or equal
  sig2.mat=matrix(sig2,nloc,ncomm-1,byrow=T)
  prec=nge[,-ncomm]+(1/sig2.mat)
  var1=1/prec
  pmedia=soma+(medias/sig2.mat)
  tmp=rnorm(nloc*(ncomm-1),mean=var1*pmedia,sd=sqrt(var1))
  psi=matrix(tmp,nloc,ncomm-1)
  
  list(psi=psi,nlk=nlk,nks=nks,soma=soma,sig2=sig2,betas=betas)
}
#---------------------------------------------
get.betas=function(psi,nloc,mu,sig2,ncomm,td=td,dtd=dtd,invT=invT,nparam=nparam){
  betas=matrix(NA,nparam,ncomm-1)
  for (i in 1:(ncomm-1)){
    prec=(1/sig2[i])*dtd+invT
    var1=solve(prec)
    pmedia=(1/sig2[i])*td%*%(psi[,i]-mu)
    betas[,i]=rmvnorm(1,var1%*%pmedia,var1)
  }
  betas
}
#---------------------------------------------
get.sig2=function(psi,nloc,mu,ncomm,dmat,betas,pot.sig2){
  p1=-(nloc/2)*log(pot.sig2)
  p2=-(1/(2*pot.sig2))
  media=mu+dmat%*%betas
  err=psi-media
  sig2=rep(NA,ncomm-1)
  for (i in 1:(ncomm-1)){
    sse=t(err[,i])%*%err[,i]
    logp=p1+p2*as.numeric(sse)
    logp1=logp-max(logp)
    tmp=exp(logp1)
    p3=tmp/sum(tmp)
    # plot(p3,type='h',main=i)
    
    #multinomial draw
    tmp=rmultinom(1,size=1,prob=p3)
    ind=which(tmp==1)
    sig2[i]=pot.sig2[ind]
  }
  sig2
}
