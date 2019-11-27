LDAcov=function(y,xmat,ncomm,phi.prior,a.gamma,b.gamma,ngibbs,nburn,var.betas){
  #----------------------------------------
  #run LDA no covariates to get initial values
  
  #get functions
  setwd('U:\\GIT_models\\git_LDA_abundance')
  source('gibbs functions.R')
  source('LDA.abundance main function.R')
  sourceCpp('aux1.cpp')

  psi=phi.prior
  gamma=0.1
  res=LDA.abundance(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,psi=psi,gamma=gamma)

  nloc=nrow(y)
  nspp=ncol(y)
  array.lsk.init=res$array.lsk
  
  #determine optimal number of groups
  nlk=apply(array.lsk.init,c(1,3),sum)
  theta=nlk/apply(nlk,1,sum)
  prop=apply(theta>0.8,2,sum) #see which communities are never above 0.8
  cond=prop!=0
  ncomm=sum(cond)
  
  #re-distribute individuals within array.lsk.init that are in eliminated communities
  array.lsk=array.lsk.init[,,cond]
  for (i in 1:nloc){
    for (j in 1:nspp){
      tmp=array.lsk.init[i,j,!cond]
      n=sum(tmp)
      if (n>0){
        z=rmultinom(1,size=n,prob=rep(1/ncomm,ncomm))
        array.lsk[i,j,]=array.lsk[i,j,]+z
      }
    }
  }
  # unique(y-apply(array.lsk,c(1,2),sum))
  
  #----------------------------------------
  #LDA with covariates
  
  #get functions
  setwd('U:\\GIT_models\\git_LDA_MS')
  source('LDA cov main function.R')
  source('LDA cov aux functions.R')
  sourceCpp('LDA_cov_aux1_cpp.cpp')
  
  res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                    phi.prior=phi.prior,array.lsk.init=array.lsk,
                    a.gamma=a.gamma,b.gamma=b.gamma,var.betas=var.betas)
  res$ncomm=ncomm
  res
}

