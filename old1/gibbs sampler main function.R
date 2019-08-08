lda.abundance.regression=function(dat,ncomm,phi.prior,ngibbs,mu,nburn,dmat){

  #useful stuff
  y=data.matrix(dat)
  nloc=nrow(y)
  nspp=ncol(y)
  nparam=ncol(dmat)
  hi=0.999999
  lo=0.000001
  td=t(dmat)
  dtd=td%*%dmat
  invT=diag(1/100,nparam)
  pot.sig2=seq(from=0.1,to=3.7,by=0.1)
  
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(1/nspp,ncomm,nspp)
  psi=matrix(0,nloc,ncomm-1)
  betas=matrix(0,nparam,ncomm-1)
  sig2=rep(3.72,ncomm-1)
  
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  psi.out=matrix(NA,ngibbs,nloc*(ncomm-1))
  beta.out=matrix(NA,ngibbs,nparam*(ncomm-1))
  sig2.out=matrix(NA,ngibbs,ncomm-1)
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
    tmp=get.psi(soma=soma,nlk=nlk,sig2=sig2,mu=mu,ncomm=ncomm,theta=theta,nks=nks,change1=change1,
                nburn=nburn,i=i,dmat=dmat,betas=betas,psi=psi)
    psi=tmp$psi
    nlk=tmp$nlk
    nks=tmp$nks
    soma=tmp$soma
    betas=tmp$betas
    sig2=tmp$sig2
    v=pnorm(psi) #calculate implied theta
    theta=convertVtoTheta(vmat=cbind(v,1),prod=rep(1,nloc))
    
    #sample phi
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    
    #sample betas and sig2
    betas=matrix(0,nparam,ncomm-1)
    if (i>(nburn/2)){  #start by fitting regular LDA and only estimate betas and sig2 afterwards
      betas=get.betas(psi=psi,nloc=nloc,mu=mu,sig2=sig2,
                      ncomm=ncomm,td=td,dtd=dtd,invT=invT,nparam=nparam)
      sig2=get.sig2(psi=psi,nloc=nloc,mu=mu,ncomm=ncomm,dmat=dmat,betas=betas,pot.sig2=pot.sig2)
    }

    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #store results  
    llk[i]=sum(y*log(prob))
    theta.out[i,]=theta
    phi.out[i,]=phi
    psi.out[i,]=psi
    beta.out[i,]=betas
    sig2.out[i,]=sig2
  }
  
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,],
       psi=psi.out[seq1,],
       betas=beta.out[seq1,],
       sig2=sig2.out[seq1,])
}
