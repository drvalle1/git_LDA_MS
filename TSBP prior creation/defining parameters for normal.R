gamma=0.1
media=1/(1+gamma); media
alpha=qnorm(media); pnorm(alpha); alpha

nsim=10000
v=rbeta(nsim,1,gamma)
mean(v==1)
v[v>0.99999999]=0.99999999
var(qnorm(v))
