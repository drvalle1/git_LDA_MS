gibbs.LDA.cov=function(ncomm,ngibbs,y,xmat,phi.prior,array.lsk.init,
                       var.betas,phi.init,estimate.phi){
  #basic settings
  nparam=ncol(xmat)
  nloc=nrow(y)
  nspp=ncol(y)
  ntot=apply(y,1,sum)
  
  #get phi by eliminating superfluous groups
  ncomm.init=ncol(phi.init)/nspp
  tmp=matrix(1:(ncomm.init*nspp),ncomm.init,nspp)
  seq1=1:ncomm
  ind1=tmp[-seq1,] #indicators for superfluous groups
  phi.mat=phi.init[,-ind1]
  phi.nrow=nrow(phi.mat)
  phi=matrix(phi.mat[phi.nrow,],ncomm,nspp)

  #get theta
  nlk=apply(array.lsk.init,c(1,3),sum)
  theta1=nlk/apply(nlk,1,sum)
  
  #re-distribute individuals within array.lsk.init that are in eliminated communities
  seq1=1:ncomm
  array.lsk=array.lsk.init[,,seq1]
  for (i in 1:nloc){
    for (j in 1:nspp){
      tmp=array.lsk.init[i,j,-seq1]
      n=sum(tmp)
      if (n>0){
        prob=theta1[i,seq1]*phi[seq1,j]
        prob=prob/sum(prob)
        z=rmultinom(1,size=n,prob=prob)
        array.lsk[i,j,]=array.lsk[i,j,]+z
      }
    }
  }

  #initial values
  nlk=apply(array.lsk,c(1,3),sum)
  betas=matrix(0,nparam,ncomm)
  options(warn=-1) #sometimes I get "glm.fit: fitted rates numerically 0 occurred" here
  for (i in 1:ncomm){
    dat.tmp=cbind(nlk[,i],xmat[,-1])
    colnames(dat.tmp)=rep('',ncol(dat.tmp)) #this is important otherwise next line breaks when we only have a single covariate
    colnames(dat.tmp)[1]='y'
    dat.tmp1=as.data.frame(dat.tmp)
    res=try(glm.nb(y ~ ., data = dat.tmp1),silent=T)
    
    #if we run into an error using NB regression, use Poisson reg
    ind=grep('Error',res)
    if (length(ind)>0) res=glm(y~.,data=dat.tmp1,family='poisson')
    
    betas[,i]=res$coef
  }
  options(warn=2)
  nks=t(apply(array.lsk,2:3,sum))
  nk=rowSums(nks)
  # phi=nks/apply(nks,1,sum); apply(phi,1,sum)
  NBN=10
  
  #to store outcomes from gibbs sampler
  phi.out=matrix(NA,ngibbs,nspp*ncomm)
  nlk.out=matrix(NA,ngibbs,nloc*ncomm)
  llk.out=rep(NA,ngibbs)
  fmodel.out=matrix(NA,ngibbs,1)
  betas.out=matrix(NA,ngibbs,nparam*ncomm)
  NBN.out=matrix(NA,ngibbs,1)
  
  #useful stuff for slice sampler algorithm
  w.betas=1
  w.NBN=10
  MaxIter=100 #to avoid overly long slice samplers
  
  #to avoid numerical issues when calculating log(p) or log(1-p)
  LoThresh=0.00000001
  UpThresh=1-LoThresh
  
  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)

    #sample NBN
    media=exp(xmat%*%betas) #get mean
    NBN=SampleNBN(Media=media,y=nlk,NBN=NBN,w=w.NBN,MaxIter=MaxIter,LoThresh=LoThresh)

    #sample betas
    betas=SampleBetas(param=betas,y=nlk,xmat=xmat,w=w.betas,nparam=nparam,
                      ncomm=ncomm,var1=var.betas,NBN=NBN,MaxIter=MaxIter,
                      LoThresh=LoThresh)

    #sample phi
    if (estimate.phi)  phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    if (!estimate.phi){
      oo=sample(phi.nrow,size=1)
      phi=matrix(phi.mat[oo,],ncomm,nspp)
    } 

    #sample z
    media=exp(xmat%*%betas) #get mean
    NBP=NBN/(media+NBN)     #get NBP
    NBP[NBP>UpThresh]=UpThresh
    tmp = SampleArray(Arraylsk=array.lsk, nloc=nloc,nspp=nspp,ncomm=ncomm,NBN=NBN,
                      y=y,LogPhi=log(phi),LogOneMinusP=log(1-NBP),
                      runif1=runif(sum(y)),nlk=nlk)
    array.lsk=tmp$ArrayLSK
    nlk=apply(array.lsk,c(1,3),sum)
    nks=t(apply(array.lsk,2:3,sum))
    nk=rowSums(nks)
    
    #calculate NB probabilities
    media.tmp=media
    media.tmp[media.tmp<LoThresh]=LoThresh
    p1=sum(dnbinom(nlk,mu=media.tmp,size=NBN,log=T))
    
    #calculate Multinom probabilities
    phi.tmp=phi
    phi.tmp[phi.tmp<LoThresh]=LoThresh
    tmp=LogLikMultin(nloc=nloc,ncomm=ncomm,nspp=nspp,LogPhi=log(phi.tmp), Arraylsk=array.lsk)
    p2=sum(tmp)

    #get phi prior
    p3=ldirichlet(x=phi.tmp,alpha=phi.prior)
    # log(ddirichlet(phi.tmp[2,],rep(phi.prior,nspp)))
    
    #get betas prior
    var.betas1=matrix(var.betas,nparam,ncomm)
    p4=dnorm(betas,mean=0,sd=sqrt(var.betas1),log=T)
    
    #store results  
    llk.out[i]=sum(p1)+sum(p2)
    fmodel.out[i]=sum(p1)+sum(p2)+sum(p3)+sum(p4)
    phi.out[i,]=phi
    nlk.out[i,]=nlk
    betas.out[i,]=betas
    NBN.out[i]=NBN
  }

  list(llk=llk.out,phi=phi.out,nlk=nlk.out,betas=betas.out,fmodel=fmodel.out,NBN=NBN.out)  
}


