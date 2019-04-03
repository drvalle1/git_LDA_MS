LDA.MS.gibbs=function(y,xmat,ncomm,ngibbs,nburnin){
  #basic settings
  nspp=ncol(y)
  nloc=nrow(y)
  npar=ncol(xmat)
  
  #useful stuff
  hi=0.999999
  lo=0.000001
  
  #initial values of parameters
  betas=matrix(0,npar,ncomm-1)
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(1/nspp,ncomm,nspp)
  
  #MH stuff
  accept1=list(betas=matrix(0,npar,ncomm-1))
  jump1=list(betas=matrix(1,npar,ncomm-1))
  accept.output=100
  
  #gibbs details
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  betas.out=matrix(NA,ngibbs,npar*(ncomm-1))
  
  llk=rep(NA,ngibbs)
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)   
    
    #sample betas
    tmp=get.betas(xmat=xmat,y=y,
                  betas=betas,phi=phi,
                  npar=npar,ncomm=ncomm,nloc=nloc,
                  jump=jump1$betas,lo=lo)
    betas=tmp$betas
    theta=tmp$theta
    accept1$betas=accept1$betas+tmp$accept
    
    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true
    
    #sample phi
    phi=rdirichlet1(alpha=nks+1,ncomm=ncomm,nspp=nspp) 
    # phi[phi>hi]=hi; phi[phi<lo]=lo
    # phi[4,]=phi.true[4,]
    
    #adapt MH
    if (i%%accept.output==0 & i<nburnin){
      tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output) 
      jump1=tmp$jump1
      accept1=tmp$accept1
    }
    
    #calculate loglikelihood
    prob=theta%*%phi; prob[prob<lo]=lo #to avoid numerical issues
    
    #store results  
    llk[i]=sum(y*log(prob))
    theta.out[i,]=theta
    phi.out[i,]=phi
    betas.out[i,]=betas
  }
  
  seq1=nburnin:ngibbs
  #output results
  list(llk=llk[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,],
       betas=betas.out[seq1,])
}


