library(MCMCpack)
library('coda')

plot(res$llk[1:ngibbs],type='l')
nburn=900
plot(res$llk[nburn:ngibbs],type='l')
#calculate effective sample size
effectiveSize(mcmc(res$betas[nburn:ngibbs,]))
plot(res$betas[nburn:ngibbs,1],type='l')

plot(res$NBN,type='l')

compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango,col='red',lwd=2)  
}

k=res$betas[ngibbs,]
nparam=ncol(xmat)
k1=matrix(k,nparam,ncomm); round(k1,2)

ordem=rep(NA,ncomm)
for (i in 2:nparam){
  ordem[i-1]=which(k1[i,]==max(k1[i,]))
}
round(k1[,ordem],2)
image(k1[,ordem])

#look at nlk
par(mfrow=c(1,1))
tmp=matrix(res$nlk[ngibbs,],nloc,ncomm); 
boxplot(tmp)
compare1(estim=jitter(tmp[,ordem]),true=jitter(nlk.true))

#look at betas
k=matrix(res$betas[ngibbs,],nparam,ncomm)
compare1(estim=k[,ordem],true=betas.true)

#look at phi
tmp=matrix(res$phi[ngibbs,],ncomm,nspp)
tmp1=tmp[ordem,]
compare1(estim=tmp1,true=phi.true)