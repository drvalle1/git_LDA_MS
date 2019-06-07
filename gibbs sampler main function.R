lda.abundance.regression=function(dat,ncomm,phi.prior,gamma1,ngibbs,sig2,sd0,mu0,nburn){

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
  mu=rep(0,ncomm-1)
  
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  psi.out=matrix(NA,ngibbs,nloc*(ncomm-1))
  mu.out=matrix(NA,ngibbs,ncomm-1)
  llk=rep(NA,ngibbs)
  
  #run gibbs sampler
  options(warn=2)
  change1=50
  for (i in 1:ngibbs){
    print(i)   

    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
    nlk=tmp$nlk
    nks=tmp$nks
    
    #get w
    tmp=get.w(nlk=nlk,psi=psi,ncomm=ncomm,nloc=nloc)
    soma=tmp$soma

    #get psi
    tmp=get.psi(soma=soma,nlk=nlk,sig2=sig2,mu=mu,ncomm=ncomm,theta=theta,nks=nks,change1=change1,nburn=nburn,i=i)
    psi=tmp$psi
    nlk=tmp$nlk
    nks=tmp$nks
    soma=tmp$soma
    mu=tmp$mu
    v=pnorm(psi) #calculate implied theta
    theta=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
    
    #sample phi
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    
    #sample mu
    mu=get.mu(psi=psi,nloc=nloc,sd0=sd0,mu0=mu0,sig2=sig2,ncomm=ncomm)
      
    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #store results  
    llk[i]=sum(y*log(prob))
    theta.out[i,]=theta
    phi.out[i,]=phi
    psi.out[i,]=psi
    mu.out[i,]=mu
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,],
       psi=psi.out[seq1,],
       mu=mu.out[seq1,])
}
