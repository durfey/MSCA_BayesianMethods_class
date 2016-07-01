##################################
# Week 4 Workshop - Notes & Code #
##################################
opar<-par(no.readonly = TRUE)


# parameters for group means
prior.mu<-0
prior.sigma<-1
sigma<-2
nGroups<-3

# simulate group means
set.seed(7)
(mus<-sort(rnorm(nGroups,mean=prior.mu,sd=prior.sigma)))
mean(mus)

# simulate data
set.seed(11)
data<-matrix(NA,500,nGroups)
for (i in 1:nGroups){
  data[,i]<-rnorm(500,mean=mus[i],sd=sigma)
}

# estimate means by averaging
(estimates1<-sort(apply(data,2,mean)))

# calc MLE
(mse1<-sum(apply(data,1,function(z) sum((z-estimates1)^2)))/length(data))

# Estimate means using fixed effects model.
ANOVA.data<-cbind(data[,1],rep(1,500))
for (i in 2:nGroups){
  ANOVA.data<-rbind(ANOVA.data,cbind(data[,i],rep(i,500)))
}
colnames(ANOVA.data)<-c("y","group")
ANOVA.data<-as.data.frame(ANOVA.data)
ANOVA.data$group<-as.factor(ANOVA.data$group)
lmod<-lm(y~.,ANOVA.data)
summary(lmod)
sort(c(lmod$coefficients[1],lmod$coefficients[-1]+lmod$coefficients[1]))



# Estimate means using random effects, which is equivalent to Bayesian estimation.
library(lme4)
lmod.ran<-lmer(y~1+(1|group),data=ANOVA.data)
(estimates2<-as.numeric(names(table(fitted(lmod.ran)))))
(mse2<-sum(apply(data,1,function(z) sum((z-estimates2)^2)))/length(data))

# show shrinkage
rbind(estimates1,estimates2)
# notice that the values from estimates2 is slightly closer to zero

c(mse1=mse1,mse2=mse2)
mean(mus)


##################
# Exchangeability and sampling

# data on 7 out of 8 states
y.7<-c(5.8,6.6,7.8,5.6,7.0,7.1,5.4)
summary(y.7)

# things we know
## 1. the 8 states are the mountain states (Arizona, Colorado, Idaho, Montana, Nevada, New Mexico, Utah, Wyoming)
## 2. Utah will have rate much lower than the others
## 3. Nevada will have rate much higher than the others

# We see that the observed rates are fairly close together, so we think the missing state is either Utah or Nevada

# If we learn that the missing state is Nevada, we can't treat the data as exchangeable since we have info that distinguishes one point (Nevada) from the others




########################
# Binomial Model with hyperparameters


##################
# High concentration, uncertain hyperparameter

# Let κ=100κ=100 be fixed and Aω=Bω=2Aω=Bω=2. These parameters correspond to density function
Omega<-Theta<-seq( 0 , 1 , length=101 )
plot(Omega,dbeta(Omega,2,2))


# Define parameters of the hyperprior.
A_omega<-2
B_omega<-2
K<-100

# The Join Prior is:
# p(θ,ω) = p(θ|ω)p(ω) = dbeta(θ|ω(100−2)+1,(1−ω)(100−2)+1)*dbeta(ω|2,2)


# Calculate and visualize joint prior and its marginals.
jointPrior<-function(theta,omega,A_omega,B_omega,K){
  res<-dbeta(omega,A_omega,B_omega)*dbeta(theta,omega*(K-2)+1,(1-omega)*(K-2)+1)
  res
}

dens<-expand.grid(Omega,Theta)
colnames(dens)<-c("Omega","Theta")
dens$Prior<-apply(dens,1,function(z) jointPrior(z[1],z[2],A_omega,B_omega,K))
Prior.theta.omega<-matrix(dens$Prior,101,101)
Prior.theta.omega<-Prior.theta.omega/sum(Prior.theta.omega) #Joint prior
Prior.omega.marginal<-apply(Prior.theta.omega,2,sum)
Prior.omega.marginal<-Prior.omega.marginal/sum(Prior.omega.marginal)*100 #Omega marginal prior
matplot(Omega,cbind(Prior.omega.marginal,dbeta(Omega,A_omega,B_omega)),type="l",ylab="Marginal p(omega)")

Prior.theta.marginal<-apply(Prior.theta.omega,1,sum)
Prior.theta.marginal<-Prior.theta.marginal/sum(Prior.theta.marginal)*100 #Theta marginal prior
plot(Theta,Prior.theta.marginal,type="l",ylab="Marginal p(theta)")

# perspective plot
persp(Theta,Omega,Prior.theta.omega,d=1,theta=-25,phi=20,main="Joint Prior Distribution")

contour(x=Omega,y=Theta,z=Prior.theta.omega,ylab="omega",xlab="theta",main="Joint Prior Distribution")


par(mfrow=c(3,1))
Prior.theta.omega.25<-jointPrior(Theta,0.25,A_omega,B_omega,K)
Prior.theta.omega.25<-Prior.theta.omega.25/sum(Prior.theta.omega.25)*100
plot(Theta,Prior.theta.omega.25,type="l",ylab="p(theta|omega=0.25)",main="Marginal prior for Theta")
Prior.theta.omega.5<-jointPrior(Theta,0.5,A_omega,B_omega,K)
Prior.theta.omega.5<-Prior.theta.omega.5/sum(Prior.theta.omega.5)*100
plot(Theta,jointPrior(Theta,0.5,A_omega,B_omega,K),type="l",ylab="p(theta|omega=0.5)")
Prior.theta.omega.75<-jointPrior(Theta,0.75,A_omega,B_omega,K)
Prior.theta.omega.75<-Prior.theta.omega.75/sum(Prior.theta.omega.75)*100
plot(Theta,jointPrior(Theta,0.75,A_omega,B_omega,K),type="l",ylab="p(theta|omega=0.75)")


for(i in seq(0.1,0.9,0.1)){
  # par(mfrow=c(3,1))
  Prior.theta.omega.i<-jointPrior(Theta,i,A_omega,B_omega,K)
  Prior.theta.omega.i<-Prior.theta.omega.i/sum(Prior.theta.omega.i)*100
  plot(Theta,Prior.theta.omega.i,type="l",ylab="p(theta|omega=0.25)",main="Marginal prior for Theta")
}
par(opar)



# likelihood function
likeli<-function(theta,s,k){
  theta^k*(1-theta)^(s-k)
}

# likelihood is based on the data of 9 successes out of 12 Bernoulli trials
likelihood<-likeli(Theta,12,9)
plot(Theta,likelihood,type="l",ylab="p(y|theta)",main="Likelihood")

Posterior <- apply(Prior.theta.omega,2,function(z)z*likelihood)
Posterior <- Posterior/sum(Posterior) # make posterior probs sum to 1
head(Posterior)

persp(Theta,Omega,Posterior,d=1,theta=-25,phi=20,main="Joint Posterior Distribution")
contour(x=Omega,y=Theta,z=Posterior,ylab="omega",xlab="theta",main="Joint Posterior Distribution")

# Calculate and plot posterior marginal distributions p(θ|y), p(ω|y)p(θ|y), p(ω|y) by adding the joint posterior matrix by row or by column.
Posterior.omega.marginal<-apply(Posterior,2,sum)
Posterior.omega.marginal<-Posterior.omega.marginal/sum(Posterior.omega.marginal)*100
Posterior.theta.marginal<-apply(Posterior,1,sum)
Posterior.theta.marginal<-Posterior.theta.marginal/sum(Posterior.theta.marginal)*100

plot(Omega,Post.omega.marginal,type="l")
plot(Theta,Post.theta.marginal,type="l")

# Show dependence of θ on ω.
par(mfrow=c(3,1))
#Omega=0.25
Post.theta.omega.25<-Posterior[,match(0.25,Omega)]
Post.theta.omega.25<-Post.theta.omega.25/sum(Post.theta.omega.25)*100
plot(Theta,Post.theta.omega.25,type="l",ylab="p(theta|omega=0.25,y)",main="Marginal posterior for Theta")
#Omega=0.5
Post.theta.omega.5<-Posterior[,match(0.5,Omega)]
Post.theta.omega.5<-Post.theta.omega.5/sum(Post.theta.omega.5)*100
plot(Theta,Post.theta.omega.5,type="l",ylab="p(theta|omega=0.25,y)",main="Marginal posterior for Theta")
#Omega=0.75
Post.theta.omega.75<-Posterior[,match(0.75,Omega)]
Post.theta.omega.75<-Post.theta.omega.75/sum(Post.theta.omega.75)*100
plot(Theta,Post.theta.omega.75,type="l",ylab="p(theta|omega=0.25,y)",main="Marginal posterior for Theta")
par(opar)


# Compare marginal priors and posteriors for both parameters.
matplot(Theta,cbind(Prior.theta.marginal,Posterior.theta.marginal),type="l",ylab="Theta Prior & Posterior")
matplot(Omega,cbind(Prior.omega.marginal,Posterior.omega.marginal),type="l",ylab="Omega Prior & Posterior")
matplot(Theta,cbind(Prior.theta.omega.25,Post.theta.omega.25),type="l",ylab="Conditional Prior and Posterior, omega=0.25")
matplot(Theta,cbind(Prior.theta.omega.5,Post.theta.omega.5),type="l",ylab="Conditional Prior and Posterior, omega=0.5")
matplot(Theta,cbind(Prior.theta.omega.75,Post.theta.omega.75),type="l",ylab="Conditional Prior and Posterior, omega=0.75")

# QUESTION: The conditional posteriors do not differ a lot from the conditional priors. What is the reason?
# ANSWER: Because we have a high concentration (K)?...




#################
# Low concentration, more certain hyperparameter