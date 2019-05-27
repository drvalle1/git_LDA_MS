lda.abundance.regression=function(dat,ncomm,phi.prior,gamma1,ngibbs,sd0,mu0,nburn){

  #useful stuff
  y=data.matrix(dat)
  nloc=nrow(y)
  nspp=ncol(y)
  hi=0.999999
  lo=0.000001
  
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(1/nspp,ncomm,nspp)
  psi=matrix(0,nloc,ncomm-1)
  
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  psi.out=matrix(NA,ngibbs,nloc*(ncomm-1))
  llk=rep(NA,ngibbs)
  
  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)   

    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
    nlk=tmp$nlk
    nks=tmp$nks
    
    #get w
    tmp=get.w(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)
    soma=tmp$soma
    nge=tmp$nge
    
    #get psi
    psi=get.psi(soma=soma,nge=nge,sd0=sd0,mu0=mu0,ncomm=ncomm)
    v=pnorm(psi) #calculate implied theta
    theta=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
    
    #sample phi
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    
    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #store results  
    llk[i]=sum(y*log(prob))
    theta.out[i,]=theta
    phi.out[i,]=phi
    psi.out[i,]=psi
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,],
       psi=psi.out[seq1,])
}
