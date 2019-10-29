LDAcov=function(y,xmat,ncomm,phi.prior,a.gamma,b.gamma,ngibbs,nburn,var.betas){
  #----------------------------------------
  #run LDA no covariates to get initial values
  
  
  #get functions
  setwd('U:\\GIT_models\\LdaPoisson_nocov')
  source('LdaPoisson_nocov main function.R')
  source('LdaPoisson_nocov aux functions.R')
  sourceCpp('LdaPoisson_nocov_aux_cpp.cpp')
  
  res=gibbs.LDA.nocov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,
                      phi.prior=phi.prior,a.gamma=a.gamma,b.gamma=b.gamma)
  
  nloc=nrow(y)
  nspp=ncol(y)
  array.lsk.init=array(res$array.lsk[ngibbs,],dim=c(nloc,nspp,ncomm))
  #----------------------------------------
  #LDA with covariates
  
  #get functions
  setwd('U:\\GIT_models\\git_LDA_MS')
  source('LDA cov main function.R')
  source('LDA cov aux functions.R')
  sourceCpp('LDA_cov_aux1_cpp.cpp')
  
  res=gibbs.LDA.cov(ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,y=y,xmat=xmat,
                    phi.prior=phi.prior,array.lsk.init=array.lsk.init,
                    a.gamma=a.gamma,b.gamma=b.gamma,var.betas=var.betas)
  res  
}

