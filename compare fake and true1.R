table(gamma.out)

boxplot(theta)

ind=1:5
ind=c(2,3,4,1,5)
theta1=theta[,ind]
plot(NA,NA,xlim=c(1,nloc),ylim=c(0,1))
for (i in 1:ncomm) lines(1:nloc,theta1[,i],col=i)

plot(phi.true,phi[ind,])