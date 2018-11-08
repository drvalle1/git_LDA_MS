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
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1){   #accept for M, M-H

  nz           <- length(x0)  #no. to accept

  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  x0[keep] <- x1[keep]           
  x0
}
#-------------------------------
#this function generates truncated normal random variables
tnorm <- function(n,lo,hi,mu,sig){   
  #normal truncated lo and hi
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo
  z[z == Inf]   <- hi
  z
}
#-----------------------------------------
#this function generates vmat, which is then used to generate the theta matrix
get.delta.theta=function(nlk,gamma,ncomm,nloc,delta,sig2,mu,jump){
  ngreater1=ngreater(nlk,nloc,ncomm)
  delta.new=delta.old=delta
  tmp=tnorm(nloc*(ncomm-1),lo=-5,hi=5,mu=delta.old,sig=jump) #need to adjust acceptance prob
  proposed=matrix(tmp,nloc,ncomm-1)
  fix.TN=fix.MH(lo=-5,hi=5,old1=delta.old,new1=proposed,jump=jump)
  
  #pre-calculate stuff
  calc.old=log(1+exp(-delta.old))
  calc.new=log(1+exp(-proposed))
  calc.old1=(1/(2*sig2))*((delta.old-mu)^2)
  calc.new1=(1/(2*sig2))*((proposed-mu)^2)
  for (i in 1:(ncomm-1)){
    delta.new=delta.old
    delta.new[,i]=proposed[,i]
  
    #calculate probabilities
    pold=-nlk[,i]*calc.old[,i]-ngreater1[,i+1]*(delta.old[,i]+calc.old[,i])-calc.old1[,i]
    pnew=-nlk[,i]*calc.new[,i]-ngreater1[,i+1]*(delta.new[,i]+calc.new[,i])-calc.new1[,i]
    
    #accept MH piece
    delta.old[,i]=acceptMH(p0=pold,p1=pnew+fix.TN[,i],x0=delta.old[,i],x1=delta.new[,i])
  }
  prob=1/(1+exp(-delta.old))
  vmat=cbind(prob,1)
  theta=convertVtoTheta(vmat,rep(1,nloc))
  list(theta=theta,delta=delta.old,accept=delta!=delta.old)  
}

#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
