loc=10
species=5

#assume that this are the results after having removed the i-th individual
int.likel.multin=function(array.lsk,prior.phi,nspp,comm.target){
  nlsk=array.lsk[,,comm.target]
  nlk=rowSums(nlsk)
  nsk=colSums(nlsk)
  nk=sum(nlsk)
  
  p1=sum(lfactorial(nlk))
  p2=sum(lfactorial(nlsk))
  p3=lgamma(nspp*prior.phi)
  p4=nspp*lgamma(prior.phi)
  p5=sum(lgamma(nsk+prior.phi))
  p6=lgamma(nk+nspp*prior.phi)
  p1-p2+p3-p4+p5-p6
}

int.likel.poisson=function(array.lsk,a.gamma,b.gamma,comm.target,xmat,betas){
  nlsk=array.lsk[,,comm.target]
  nlk=rowSums(nlsk)

  lmedia=xmat%*%betas[,comm.target]
  
  p1=sum(lmedia)
  p2=sum(lfactorial(nlk))
  p3=a.gamma*log(b.gamma)
  p4=lgamma(a.gamma)
  p5=lgamma(a.gamma+sum(nlk))
  tmp=sum(exp(lmedia))+b.gamma
  p6=(a.gamma+sum(nlk))*log(tmp)

  p1-p2+p3-p4+p5-p6
}

#calculate using full equations
prob=rep(0,ncomm)
for (ksel in 1:ncomm){ 
  #adding one individual in this group
  array.lsk.new=array.lsk
  array.lsk.new[loc,species,ksel]=array.lsk.new[loc,species,ksel]+1

  #get multinomial and poisson integrated likelihood
  lmultin=lpoisson=rep(NA,ncomm)
  for (j in 1:ncomm){
    lmultin[j]=int.likel.multin(array.lsk=array.lsk.new,prior.phi=prior.phi,nspp=nspp,comm.target=j)
    lpoisson[j]=int.likel.poisson(array.lsk=array.lsk.new,a.gamma=0.1,b.gamma=0.1,comm.target=j,
                                  xmat=xmat,betas=betas)
  }
  prob[ksel]=sum(lmultin)+sum(lpoisson)
}
prob=prob-max(prob)
prob=exp(prob)
prob=prob/sum(prob)
prob.full=prob

#calculate using subset of expressions
prob=rep(0,ncomm)
for (ksel in 1:ncomm){
  lmedia=xmat[loc,]%*%betas[,ksel]
  p1=lmedia
  media=exp(lmedia)
  p2=log(b.gamma+sum(media))
  p3=log(sum(array.lsk[,species,ksel])+prior.phi)
  nk=sum(array.lsk[,,ksel])
  p4=log(nk+nspp*prior.phi)
  p5=log(a.gamma+nk)
  p6=log(array.lsk[loc,species,ksel]+1)
  prob[ksel]=p1-p2+p3-p4+p5-p6
}
prob=prob-max(prob)
prob=exp(prob)
prob=prob/sum(prob)
prob.subset=prob
