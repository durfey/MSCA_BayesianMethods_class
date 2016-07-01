datapath<-"./DBDA2Eprograms"

source(paste(dataPath,"DBDA2E-utilities.R",sep="/"))  # Load definitions of graphics functions etc.
source('./DBDA2Eprograms/DBDA2E-utilities.R')
source('./DBDA2Eprograms/BernBeta.R')




# t<-omega<-0.5 # prior mode
# k<-kappa<-500 # concentration
# a<-b<-alpha<-beta<-250 # from belief that coin is fair (symmetric)
# s<-17
# 
# a<-
# b<-(1-t)*(k-2)+1



# Specify the prior:
t = 0.5             # Specify the prior MODE.
n = 500               # Specify the effective prior sample size.
a = t*(n-2) + 1      # Convert to beta shape parameter a.
b = (1-t)*(n-2) + 1  # Convert to beta shape parameter b.

(Prior = c(a,b))  




# Specify the data:
N = 20                         # The total number of flips.
z = 17                         # The number of heads.
Data = c(rep(0,N-z),rep(1,z))  # Convert N and z into vector of 0's and 1's.


# Use the same function for calculating binomial likelihood as in the previous workshop.
likeli<-function(par,data){
  sdata<-sum(data)
  ldata<-length(data)
  return(par^sdata*(1-par)^(ldata-sdata))
}


Theta = seq( 0 , 1 , length=1001 )  # Make grid of Theta.
pTheta = dbeta(Theta ,shape1=Prior[1],shape2=Prior[2]) # Beta prior distribution.
plot(Theta,dbeta(Theta ,shape1=Prior[1],shape2=Prior[2]),ylab="Prior")


likelihood<-likeli(Theta,Data)
plot(Theta,likelihood)


post<-pTheta*likelihood
post<-post/sum(post)*1000


plot (Theta,post)


(mode<-Theta[which.max(post)])


#openGraph(width=5,height=7)
posterior = BernBeta( priorBetaAB=Prior, Data=Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )







Prior<-c(a=1,b=1)
Data<-c(s=100,k=58)
Posterior<-c(Prior["a"]+Data["k"],Prior["b"]+Data["s"]-Data["k"])
(HDI<-c(qbeta(.025,Posterior["a"],Posterior["b"]),
        qbeta(.975,Posterior["a"],Posterior["b"])))






kappa<-1000
mu<-.5
alpha<-mu*kappa
beta<-(1-mu)*kappa
data<-c(rep(1,9),rep(0,1))
data.k<-sum(data)
data.s<-length(data)
alpha.post<-alpha+data.k
beta.post<-data.s-data.k+beta
mean.posterior<-alpha.post/(alpha.post+beta.post)
likelihood.postTheta<-mean.posterior^1*(1-mean.posterior)^0