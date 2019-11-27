gibbs.LDA.cov=function(ncomm,ngibbs,nburn,y,xmat,phi.prior,array.lsk.init,
                       a.gamma,b.gamma,var.betas){
  #basic settings
  nparam=ncol(xmat)
  nloc=nrow(y)
  nspp=ncol(y)
  
  #initial values
  array.lsk=array.lsk.init
  nlk=apply(array.lsk,c(1,3),sum)
  betas=matrix(0,nparam,ncomm)
  lambda=rep(0,ncomm)
  options(warn=-1) #sometimes I get "glm.fit: fitted rates numerically 0 occurred" here
  for (i in 1:ncomm){
    dat.tmp=cbind(nlk[,i],xmat)
    colnames(dat.tmp)=c('y',paste0('cov',1:ncol(xmat)))
    dat.tmp1=as.data.frame(dat.tmp)
    res=glm(y~.,data=dat.tmp1,family='poisson')
    lambda[i]=exp(res$coef[1])
    betas[,i]=res$coef[-1]
  }
  options(warn=2)
  nks=t(apply(array.lsk,2:3,sum))
  phi=nks/apply(nks,1,sum); apply(phi,1,sum)

  #to store outcomes from gibbs sampler
  lambda.out=matrix(NA,ngibbs,ncomm)
  phi.out=matrix(NA,ngibbs,nspp*ncomm)
  nlk.out=matrix(NA,ngibbs,nloc*ncomm)
  llk.out=rep(NA,ngibbs)
  llk.ind.out=matrix(NA,ngibbs,nloc)
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
    nlk=apply(array.lsk,c(1,3),sum)
    nks=t(apply(array.lsk,2:3,sum))

    #sample betas
    tmp=sample.betas(nlk=nlk,xmat=xmat,betas=betas,
                     ncomm=ncomm,nparam=nparam,jump=jump1$betas,
                     a.gamma=a.gamma,b.gamma=b.gamma,var.betas=var.betas)
    betas=tmp$betas
    accept1$betas=accept1$betas+tmp$accept
    
    #sample phi
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)

    #sample lambda
    lambda=get.lambda(nlk=nlk,a.gamma=a.gamma,b.gamma=b.gamma,
                      xmat=xmat,betas=betas,ncomm=ncomm)

    #adaptive MH
    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #calculate Poisson probabilities
    media=matrix(lambda,nloc,ncomm,byrow=T)*exp(xmat%*%betas)
    p1=dpois(nlk,media,log=T)
    
    #calculate loglik
    phi.tmp=phi; phi.tmp[phi.tmp<0.00001]=0.00001
    p2=LogLikMultin(nloc=nloc,ncomm=ncomm,nspp=nspp,phi=phi.tmp,Arraylsk=array.lsk)    
    
    #get phi prior
    p3=ldirichlet(x=phi.tmp,alpha=phi.prior)
    # log(ddirichlet(phi[2,],rep(phi.prior,nspp)))
    
    #get betas prior
    p4=dnorm(betas,mean=0,sd=sqrt(var.betas),log=T)
    
    #get lambda prior
    p5=dgamma(lambda,a.gamma,b.gamma,log=T)
    
    #store results  
    llk.out[i]=sum(p1)+sum(p2)+sum(p3)+sum(p4)+sum(p5)
    llk.ind.out[i,]=p2 
    phi.out[i,]=phi
    lambda.out[i,]=lambda
    nlk.out[i,]=nlk
    betas.out[i,]=betas
  }
  
  list(llk=llk.out,phi=phi.out,lambda=lambda.out,nlk=nlk.out,betas=betas.out,llk.ind.out=llk.ind.out)  
}


