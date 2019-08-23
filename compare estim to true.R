plot(llk.out[100:ngibbs],type='l')

compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango)  
}

ordem=c(4,2,3,1,5)

#look at lambda
compare1(estim=lambda.out[ngibbs,ordem],true=lambda.true)

#look at nlk
tmp=matrix(nlk.out[ngibbs,],nloc,ncomm)
compare1(estim=tmp[,ordem],true=nlk.true)

#look at betas
k=matrix(betas.out[ngibbs,],nparam,ncomm)
compare1(estim=k[,ordem],true=betas.true)

#look at phi
tmp=matrix(phi.out[ngibbs,],ncomm,nspp)
tmp1=tmp[ordem,]
compare1(estim=tmp1,true=phi.true)

#look at convergence
par(mfrow=c(4,3),mar=rep(1,4))
for (i in 1:12){
  plot(betas.out[,i],type='l')
}

par(mfrow=c(3,2),mar=rep(1,4))
for (i in 1:ncomm){
  plot(lambda.out[,i],type='l')
}
