gibbs.LDA.cov=function(ncomm,ngibbs,nburn,y,xmat,phi.prior,array.lsk.init,
                       var.betas){
  #basic settings
  nparam=ncol(xmat)
  nloc=nrow(y)
  nspp=ncol(y)
  ntot=apply(y,1,sum)
  
  #initial values
  array.lsk=array.lsk.init
  nlk=apply(array.lsk,c(1,3),sum)
  betas=matrix(0,nparam,ncomm)
  options(warn=-1) #sometimes I get "glm.fit: fitted rates numerically 0 occurred" here
  for (i in 1:ncomm){
    dat.tmp=cbind(nlk[,i],xmat[,-1])
    colnames(dat.tmp)[1]='y'
    dat.tmp1=as.data.frame(dat.tmp)
    res=glm(y~.,data=dat.tmp1,family='poisson')
    betas[,i]=res$coef
  }
  options(warn=2)
  nks=t(apply(array.lsk,2:3,sum))
  phi=nks/apply(nks,1,sum); apply(phi,1,sum)

  #to store outcomes from gibbs sampler
  phi.out=matrix(NA,ngibbs,nspp*ncomm)
  nlk.out=matrix(NA,ngibbs,nloc*ncomm)
  llk.out=rep(NA,ngibbs)
  fmodel.out=matrix(NA,ngibbs,1)
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

    #get summaries of array.lsk
    nk=colSums(nlk)

    #sample z
    tmp = SampleArray(Arraylsk=array.lsk, nloc=nloc,nspp=nspp,ncomm=ncomm,
                      y=y,lpmedia=lpmedia,runif1=runif(sum(y)),
                      nk=nk,nks=nks, PriorPhi=phi.prior)
    array.lsk=tmp$ArrayLSK
    # array.lsk=array.lsk.true
    nlk=apply(array.lsk,c(1,3),sum)
    nks=t(apply(array.lsk,2:3,sum))

    #sample betas
    tmp=sample.betas(y=y,xmat=xmat,betas=betas,
                     ncomm=ncomm,nparam=nparam,jump=jump1$betas,
                     var.betas=var.betas,phi=phi,ntot=ntot)
    betas=tmp$betas
    accept1$betas=accept1$betas+tmp$accept
    # betas=rbind(lambda.true,betas.true)
    
    #sample phi
    phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    # phi=phi.true
    
    #adaptive MH
    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #calculate Poisson probabilities
    media=exp(xmat%*%betas)
    p1=dpois(ntot,apply(media,1,sum),log=T)
    
    #calculate loglik
    theta=media/apply(media,1,sum)
    prob=theta%*%phi
    p2=y*log(prob)
    
    #get phi prior
    phi.tmp=phi; phi.tmp[phi.tmp<0.00000000001]=0.00000000001
    p3=ldirichlet(x=phi.tmp,alpha=phi.prior)
    # log(ddirichlet(phi[2,],rep(phi.prior,nspp)))
    
    #get betas prior
    var.betas1=matrix(var.betas,nparam,ncomm)
    p4=dnorm(betas,mean=0,sd=sqrt(var.betas1),log=T)
    
    #store results  
    llk.out[i]=sum(p1)+sum(p2)
    fmodel.out[i]=sum(p1)+sum(p2)+sum(p3)+sum(p4)
    phi.out[i,]=phi
    nlk.out[i,]=nlk
    betas.out[i,]=betas
  }
  
  list(llk=llk.out,phi=phi.out,nlk=nlk.out,betas=betas.out,fmodel=fmodel.out)  
}


