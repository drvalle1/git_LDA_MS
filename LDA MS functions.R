#this function adjusts the jump parameter in MH algorithm
#with the goal of achieving an acceptance rate between 0.2 and 0.4
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<5
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
#accepts or rejects based on MH algorithm
acceptMH <- function(p0,p1,x0,x1){   #accept for M, M-H

  nz           <- length(x0)  #no. to accept

  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  x0[keep] <- x1[keep]           
  x0
}
#-------------------------------
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
#-----------------------------------------
get.theta=function(media,nloc){
  probs=1/(1+exp(-media))
  vmat=cbind(probs,1)
  convertVtoTheta(vmat,rep(1,nloc))
}
#-----------------------------------------
#this function samples the regression slope parameters
get.betas=function(xmat,y,betas,phi,npar,ncomm,jump,lo,nloc){
  betas.old=betas.new=betas.orig=betas
  media.old=xmat%*%betas.old
  
  #propose new beta values
  tmp=rnorm(npar*(ncomm-1),mean=betas.old,sd=jump)
  proposed=matrix(tmp,npar,ncomm-1)
  
  #calculate prior probabilities
  prior.old=(-1/2)*(betas.old^2)
  prior.new=(-1/2)*(proposed^2)
  
  for (i in 1:npar){
    for (j in 1:(ncomm-1)){
      betas.new=betas.old
      media.new=media.old
      
      betas.new[i,j]=proposed[i,j]
      media.new[,j]=xmat%*%betas.new[,j]
      
      theta.old=get.theta(media=media.old,nloc=nloc)
      theta.new=get.theta(media=media.new,nloc=nloc)
      prob.old=theta.old%*%phi; prob.old[prob.old<lo]=lo #to avoid numerical issues
      prob.new=theta.new%*%phi; prob.new[prob.new<lo]=lo #to avoid numerical issues
      lprob.old=log(prob.old)
      lprob.new=log(prob.new)
      pold=sum(y*lprob.old)
      pnew=sum(y*lprob.new)
      k=acceptMH(p0=pold+prior.old[i,j],
                 p1=pnew+prior.new[i,j],
                 x0=betas.old[i,j],
                 x1=betas.new[i,j])
      betas.old[i,j]=k
      media.old[,j]=xmat%*%betas.old[,j]
    }
  }
  media=xmat%*%betas.old
  theta=get.theta(media=media,nloc=nloc)
  list(betas=betas.old,
       accept=betas.old!=betas.orig,
       theta=theta)
}

