sample.betas.R=function(y,xmat,betas,ncomm,nparam,jump,var.betas,phi,ntot){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump)

  var.betas1=matrix(var.betas,nparam,ncomm)
  prior.old=dnorm(betas.orig,mean=0,sd=sqrt(var.betas1),log=T)
  prior.new=dnorm(betas.prop,mean=0,sd=sqrt(var.betas1),log=T)

  media.old=media.new=exp(xmat%*%betas.old)
  for (i in 1:nparam){
    for (j in 1:ncomm){
      betas.new=betas.old
      betas.new[i,j]=betas.prop[i,j]
      media.new[,j]=exp(xmat%*%betas.new[,j])
      soma.media.old=rowSums(media.old)
      soma.media.new=rowSums(media.new)
      theta.old=media.old/soma.media.old
      theta.new=media.new/soma.media.new

      prob.old=theta.old%*%phi
      prob.new=theta.new%*%phi
      p1.old=sum(y*log(prob.old))
      p1.new=sum(y*log(prob.new))
      p2.old=sum(dpois(ntot,soma.media.old,log=T))
      p2.new=sum(dpois(ntot,soma.media.new,log=T))

      pold=p1.old+p2.old+prior.old[i,j]
      pnew=p1.new+p2.new+prior.new[i,j]
      pthresh=exp(pnew-pold)
      if (runif(1)<pthresh){
        betas.old[i,j]=betas.new[i,j]
        media.old[,j]=media.new[,j]
      }
    }
  }

  list(betas=betas.old,accept=betas.old!=betas.orig)
}
#---------------------------------
sample.betas.arma=function(y,xmat,betas,ncomm,nparam,jump,var.betas,phi,ntot){
  betas.orig=betas.old=betas.prop=betas
  betas.prop[]=rnorm(nparam*ncomm,mean=betas.old,sd=jump)
  
  var.betas1=matrix(var.betas,nparam,ncomm)
  prior.old=dnorm(betas.orig,mean=0,sd=sqrt(var.betas1),log=T)
  prior.new=dnorm(betas.prop,mean=0,sd=sqrt(var.betas1),log=T)
  
  runif1=matrix(runif(nparam*ncomm),nparam,ncomm)
  betas=sampleBetas(y=y,xmat=xmat,betas_prop=betas.prop,
                    prior_old=prior.old, prior_new=prior.new, runif1=runif1,
                    betas=betas, phi=phi,ntot=ntot,
                    ncomm=ncomm,nparam=nparam, nloc=nloc, nspp=nspp)
  
  list(betas=betas,accept=betas!=betas.orig)
}

microbenchmark(
  sample.betas.arma(y=y,xmat=xmat,betas=betas,
                   ncomm=ncomm,nparam=nparam,jump=jump1$betas,
                   var.betas=var.betas,phi=phi,ntot=ntot),
  sample.betas.R(y=y,xmat=xmat,betas=betas,
                 ncomm=ncomm,nparam=nparam,jump=jump1$betas,
                 var.betas=var.betas,phi=phi,ntot=ntot),
  times=10
)
