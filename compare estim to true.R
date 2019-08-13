plot(llk.out[100:ngibbs],type='l')

compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango)  
}

#we need to re-order things first
compare1(estim=lambda.out[ngibbs,],true=lambda.true)
compare1(estim=nlk.out[ngibbs,],true=nlk.true)
compare1(estim=betas.out[ngibbs,],true=betas.true)

par(mfrow=c(4,3),mar=rep(1,4))
for (i in 1:12){
  plot(betas.out[,i],type='l')
}
par(mfrow=c(3,2),mar=rep(1,4))
for (i in 1:ncomm){
  plot(lambda.out[,i],type='l')
}

compare1(estim=phi.out[ngibbs,],true=phi.true)