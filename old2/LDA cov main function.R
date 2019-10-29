gibbs.LDA.cov=function(ncomm,ngibbs,nburn,y,xmat,phi.prior){
  #basic settings
  nparam=ncol(xmat)
  nloc=nrow(y)
  nspp=ncol(y)
  
  #initial values
  betas=matrix(0,nparam,ncomm)
  array.lsk=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:nloc){
    for (j in 1:nspp){
      if (y[i,j]!=0){
        array.lsk[i,j,]=rmultinom(1,size=y[i,j],prob=rep(1/ncomm,ncomm))  
      }
    }
  }

  #basic test
  # z=apply(array.lsk,1:2,sum)
  # unique(y-z)
  
  nlk=apply(array.lsk,c(1,3),sum)
  nks=t(apply(array.lsk,2:3,sum))

  phi=matrix(1/nspp,ncomm,nspp)  
  lambda=apply(nlk,2,mean)
  #get thetas
  theta=get.theta.from.lambda(lambda=lambda,ncomm=ncomm)

  #just checking
  # teste=GetLambda(LogTheta=log(theta),ncomm=ncomm)
  
  #priors
  a.gamma=b.gamma=0.1
  
  #to store outcomes from gibbs sampler
  lambda.out=matrix(NA,ngibbs,ncomm)
  phi.out=matrix(NA,ngibbs,nspp*ncomm)
  nlk.out=matrix(NA,ngibbs,nloc*ncomm)
  llk.out=rep(NA,ngibbs)
  betas.out=matrix(NA,ngibbs,nparam*ncomm)
  
  #useful stuff for MH algorithm
  accept1=list(betas=matrix(0,nparam,ncomm))
  jump1=list(betas=matrix(0.1,nparam,ncomm))
  accept.output=50
  nadapt=ngibbs/2
  
  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)

    #order communities according to size
    # if (i<nburn & i%%50==0){
    #   ind=order(lambda,decreasing=T)
    #   lambda=lambda[ind]
    #   theta=get.theta.from.lambda(lambda=lambda,ncomm=ncomm)
    #   
    #   betas=betas[,ind]
    #   phi=phi[ind,]
    #   array.lsk=array.lsk[,,ind]
    #   nlk=nlk[,ind]
    # }
    
    #get log of part 1
    lpmedia=xmat%*%betas
    pmedia=exp(lpmedia)
    pmedia.soma=colSums(pmedia)+b.gamma
    lp1=lpmedia-matrix(log(pmedia.soma),nloc,ncomm,byrow=T)
    
    #get summaries of array.lsk
    nk=colSums(nlk)
    
    #sample z
    tmp = SampleArray(Arraylsk=array.lsk, nloc=nloc,nspp=nspp,ncomm=ncomm,
                      y=y,lp1=lp1,runif1=runif(sum(y)),
                      nk=nk,nks=nks, PriorPhi=phi.prior, agamma=a.gamma)
    array.lsk=tmp$ArrayLSK
    # array.lsk=array.lsk.true
    nlk=apply(array.lsk,c(1,3),sum)
    nks=t(apply(array.lsk,2:3,sum))
    # nks=nks.true#rbind(nks.true,0,0)
    # nlk=nlk.true#cbind(nlk.true,0,0)

    #sample betas
    tmp=sample.betas(nlk=nlk,xmat=xmat,betas=betas,
                     ncomm=ncomm,nparam=nparam,jump1=jump1$betas,
                     a.gamma=a.gamma,b.gamma=b.gamma)
    betas=tmp$betas
    accept1$betas=accept1$betas+tmp$accept
    # betas=betas.true#cbind(betas.true,0,0)

    #sample thetas
    # theta=sample.theta(theta=theta,a1=a1,b1=b1,a2=a2,b2=b2,nlk=nlk,ncomm=ncomm,nloc=nloc,
    #                    xmat=xmat,betas=betas)
    # lambda=GetLambda(LogTheta=log(theta),ncomm=ncomm)
    lambda=lambda.true#c(lambda.true,0.01,0.01)
    
    #sample phi
    # phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    # phi=phi.true #rbind(phi.true,0,0)

    #adaptive MH
    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #calculate Poisson probabilities
    media=matrix(lambda,nloc,ncomm,byrow=T)*exp(xmat%*%betas)
    p1=dpois(nlk,media,log=T)
    # phi.tmp=phi; phi.tmp[phi.tmp<0.00001]=0.00001
    
    p2=LogLikMultin(nloc=nloc,ncomm=ncomm,nspp=nspp,phi=phi,Arraylsk=array.lsk)    
    
    #get phi prior
    p3=ldirichlet(x=phi,alpha=phi.prior)
    # log(ddirichlet(phi[2,],rep(phi.prior,nspp)))
    
    #get betas prior
    p4=dnorm(betas,mean=0,sd=1,log=T)
    
    #get lambda prior
    # p5=dgamma(theta[1],a1,b1,log=T)+sum(dgamma(theta[-1],a2,b2,log=T))
    
    #store results  
    llk.out[i]=sum(p1)+sum(p2)+sum(p3)+sum(p4)#+p5
    phi.out[i,]=phi
    lambda.out[i,]=lambda
    nlk.out[i,]=nlk
    betas.out[i,]=betas
  }

  list(llk=llk.out,phi=phi.out,lambda=lambda.out,nlk=nlk.out,betas=betas.out)  
}


