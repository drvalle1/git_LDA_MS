LDA.MS.gibbs=function(y,xmat,ncomm,ngibbs,nburnin,vmat.init=NULL,phi.init=NULL,phi.prior){
  #basic settings
  nspp=ncol(y)
  nloc=nrow(y)
  npar=ncol(xmat)
  nind.loc=apply(y,1,sum)
  loc.expansion=rep(1:nloc,times=nind.loc)
  
  #useful stuff
  hi=0.999999
  lo=0.000001
  tx=t(xmat)
  xtx=tx%*%xmat

  phi=matrix(1/nspp,ncomm,nspp)
  theta=matrix(1/ncomm,nloc,ncomm)
  betas=matrix(0,npar,ncomm-1)
  #initial values of parameters
  if (!is.null(vmat.init)) {
    vmat.init[vmat.init<lo]=lo
    vmat.init[vmat.init>hi]=hi
    theta=convertVtoTheta(vmat.init,rep(1,nloc))

    #calculate implied betas
    tmp=qnorm(vmat.init[,-ncomm])
    betas=solve(t(xmat)%*%xmat)%*%t(xmat)%*%tmp #these should be good starting values
  }
  if (!is.null(phi.init))   phi=phi.init

  #gibbs details
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  betas.out=matrix(NA,ngibbs,npar*(ncomm-1))
  llk=rep(NA,ngibbs)
  
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)   
    
    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true
    
    #sample phi
    # phi=rbind(phi.true,matrix(1/nspp,11,nspp))
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    # # phi[phi>hi]=hi; phi[phi<lo]=lo
    # # phi[4,]=phi.true[4,]

    #sample w
    tmp=get.w(xmat=xmat,nlk=nlk,betas=betas,ncomm=ncomm,nloc=nloc)    
    soma=tmp$soma
    n=tmp$n

    #sample betas
    tmp=get.betas(tx=tx,xmat=xmat,ncomm=ncomm,npar=npar,n=n,soma=soma,nloc=nloc)
    betas=tmp$betas
    theta=tmp$theta

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


