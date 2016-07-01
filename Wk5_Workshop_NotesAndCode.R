################################
# Week 5 Workshop Notes & Code #
################################



probabilities<-c(.1,.3,.5,.05,.05)
sum(probabilities)

init<-2

# function for making one step based on probabilities -- MY CODE ATTEMPT
onestep.mine<-function(probabilities=probabilities,position=init){
  if(init ==1)
    move.proposed <- sample(c(position+1,length(probabilities)),1,prob=c(0.5,0.5))
  else if(init == length(probabilities))
    move.proposed <- sample(c(position-1,1),1,prob=c(0.5,0.5))
  else
    move.proposed <- sample(c(position+1,position-1),1,prob=c(0.5,0.5))
  
  prob.current <- probabilities[position]
  prob.proposed <- probabilities[move.proposed]
  prob.move <- min(prob.proposed/prob.current,1)
  move.accepted <- sample(c(0,1),1,prob=c(1-prob.move,prob.move))
  if(move.accepted==1)
    position <- move.proposed
  position
}

# CODE FROM WORKSHOP
oneStep<-function(probabilities){
  numberStates<-length(probabilities)
  turn<-sample(c(-1,1),1)
  newStates<-c(init, init+turn)
  if (newStates[2]==numberStates+1) newStates[2]<-1
  if (newStates[2]==0) newStates[2]<-numberStates
  if (probabilities[newStates[2]]>probabilities[newStates[1]]) { 
    init<-newStates[2]
  } else {  
    p<-probabilities[newStates[2]]/probabilities[newStates[1]]
    init<-sample(newStates,1,prob = c(1-p,p))
  }  
  init
}  

trajectory<-rep(NA,100000)


# using my function
for (i in 1:100000){
  init<-onestep.mine(probabilities)
  trajectory[i]<-init
}
table(trajectory)
table(trajectory)/sum(table(trajectory))

# using given function
for (i in 1:100000){
  init<-oneStep(probabilities)
  trajectory[i]<-init
}
table(trajectory)
table(trajectory)/sum(table(trajectory))



# using JAGS
dataPath<-'/Users/rdurfey/R_misc/BayesianMethods/DBDA2Eprograms'
source(paste(dataPath,"DBDA2E-utilities.R",sep="/"))

myData<-read.csv(paste(dataPath,"z15N50.csv",sep="/"))
head(myData)

y<-myData$y
(Ntotal<-length(y))
(dataList<-list(y=y,Ntotal=Ntotal))

# Preparing model for JAGS
modelString=" 
model {
for (i in 1:Ntotal) {
y[i]~dbern(theta)
}
theta~dbeta(1,1)
}
"
writeLines(modelString,con="Tempmodel.txt")

# Initializing Markov Chains
MLE<-sum(y)/Ntotal
init1<-MLE
init2<-MLE*(1+.01)
init3<-MLE*(1-.01)
initsList<-function(){
  thetaInit<-sample(c(init1,init2,init3),1,replace=T)
  return(list(theta=thetaInit))
}
initsList()

# Sending information to JAGS
library(rjags)
jagsModel<-jags.model(file="TempModel.txt",data=dataList,n.chains=3,n.adapt=500)
names(jagsModel)

# Running MCMC on JAGS
update(jagsModel,n.iter=600)

codaSamples<-coda.samples(jagsModel,variable.names=c("theta"),n.iter=3334)
list.samplers(jagsModel)
head(codaSamples)
summary(codaSamples)

traceplot(codaSamples)
densplot(codaSamples)
plot(codaSamples)
autocorr.plot(codaSamples,ask=F)
effectiveSize(codaSamples)
gelman.diag(codaSamples)
gelman.plot(codaSamples)

lapply(codaSamples,mean)
sum(y)/Ntotal

(l<-min(unlist(codaSamples))-.05)
(h<-max(unlist(codaSamples))+.05)

histBreaks<-seq(l,h,by=.05)
postHist<-lapply(codaSamples,hist,breaks=histBreaks)

plot(postHist[[1]]$mids,postHist[[1]]$density,type="l",col="black",lwd=2,ylim=c(0,6),ylab="Distribution Density",xlab="Theta")
lines(postHist[[2]]$mids,postHist[[2]]$density,type="l",col="red",lwd=2)
lines(postHist[[3]]$mids,postHist[[3]]$density,type="l",col="blue",lwd=2)
lines(postHist[[3]]$mids,dbeta(postHist[[3]]$mids,1+sum(y),Ntotal-sum(y)+1),type="l",col="green",lwd=3)
legend("topright",legend=c("Chain1","Chain2","Chain3","Theoretical"),col=c("black","red","blue","green"),lwd=2)


##### section 4 #####

# Comparison of 2 Binomial 
Theta<-seq(0,1,length=1001)
plot(Theta,dbeta(Theta,1,1))

myData<-read.csv("//Users/rdurfey/R_misc/BayesianMethods/Week5/2GroupsStudy.csv")
(y=myData$y)
(s<-as.numeric(myData$s))
(Ntotal<-length(y))
(Nsubj<-length(unique(s)))
(dataList<-list(y=y,s=s,Ntotal=Ntotal,Nsubj=Nsubj))

modelString = "
  model {
for ( i in 1:Ntotal ) {
y[i] ~ dbern( theta[s[i]] )
}
for ( sIdx in 1:Nsubj ) {
theta[sIdx] ~ dbeta( 2 , 2 ) # N.B.: 2,2 prior; change as appropriate.
}
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

initsList = function() {
  thetaInit = rep(0,Nsubj)
  for ( sIdx in 1:Nsubj ) { # for each subject
    includeRows = ( s == sIdx ) # identify rows of this group
    yThisSubj = y[includeRows]  # extract data of this group
    resampledY = sample( yThisSubj , replace=TRUE ) # resample
    thetaInit[sIdx] = sum(resampledY)/length(resampledY) 
  }
  thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
  return( list( theta=thetaInit ) )
}


parameters = c( "theta")     # The parameters to be monitored
adaptSteps = 500             # Number of steps to adapt the samplers
burnInSteps = 500            # Number of steps to burn-in the chains
nChains = 4                  # nChains should be 2 or more for diagnostics 
numSavedSteps<-50000
nIter = ceiling(numSavedSteps / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )

update( jagsModel , n.iter=burnInSteps )

codaSamples = coda.samples( jagsModel , variable.names=parameters , n.iter=nIter)
head(codaSamples)

list.samplers(jagsModel)

summary(codaSamples)
plot(codaSamples)
autocorr.plot(codaSamples,ask=F)
effectiveSize(codaSamples)
gelman.diag(codaSamples)
gelman.plot(codaSamples)

matrix(unlist(lapply(codaSamples,function(z) apply(z,2,mean))),ncol=2,byrow = T)

plot(density(codaSamples[[1]][,1]),xlim=c(0,1),ylim=c(0,3))
lines(density(codaSamples[[1]][,2]))
lines(density(codaSamples[[2]][,1]),col="red")
lines(density(codaSamples[[2]][,2]),col="red")
lines(density(codaSamples[[3]][,1]),col="blue")
lines(density(codaSamples[[3]][,2]),col="blue")
lines(density(codaSamples[[4]][,1]),col="green")
lines(density(codaSamples[[4]][,2]),col="green")