################################
# Week 6 Workshop Notes & Code #
################################


##### 1. Model Comparison #####

modelString=" 
model {
for (i in 1:Ntotal) {
y[i]~dbern(theta)
}
theta~dbeta(omega[m]*(kappa-2)+1,(1-omega[m])*(kappa-2)+1)
omega[1]<-.25
omega[2]<-.75
kappa<-12
m~dcat(mPriorProb[])
mPriorProb[1]<-.5
mPriorProb[2]<-.5
}
"
writeLines(modelString,con="Tempmodel.txt")

# Create list of data corresponding to 6 successes out of 9 trials.
y<-c(rep(0,3),rep(1,6))
(Ntotal<-length(y))
(dataList<-list(y=y,Ntotal=Ntotal))

library(rjags)

jagsModel<-jags.model(file="TempModel.txt",data=dataList,n.chains=4,n.adapt=500)
names(jagsModel)

# Burn in
update(jagsModel,n.iter=600)

# Generate MCMC trajectories.
codaSamples<-coda.samples(jagsModel,variable.names=c("m"),thin=1,n.iter=5000)
list.samplers(jagsModel)
head(codaSamples)

# Analyze convergence.
summary(codaSamples)
plot(codaSamples)
autocorr.plot(codaSamples,ask=F)
effectiveSize(codaSamples)
# We see that autocorrelation converges to zero only about 5 lags or so.
# This is confirmed by ESS. Here, ESS is roughly one fifth of what was simulated

# Rerun the chains with thinning paramater equal to 5.
codaSamples<-coda.samples(jagsModel,variable.names=c("m"),thin=5,n.iter=5000)
plot(codaSamples)
autocorr.plot(codaSamples,ask=F)
effectiveSize(codaSamples)
lapply(codaSamples,length)
# Now autocorrelation function is not significant.
# Effective size is 3274.6461639, but this is out of total 4,000 of observations instead of 20,000.

# Potential scale reduction factor or shrink factor showed convergence.
gelman.diag(codaSamples)
gelman.plot(codaSamples)

# Now analyze the results.
# Look at the chain means.
(means<-lapply(codaSamples,mean))

# Find posterior probabilities of m=2m=2 for each of the 4 chains and their average.
(prob.2<-lapply(means,function(z) z-1))
mean(unlist(prob.2))
# This means that posterior probability of m=1m=1 is 0.17225 (or whatever one minus above is).

# Find how much time each chain spent in each of the state for mm.
lapply(codaSamples,function(z) sum(z==2)/length(z)) # for how much time spent in state m=2
lapply(codaSamples,function(z) sum(z==1)/length(z)) # for how much time spent in state m=1




##### 2. Application of Stan #####

source('/Users/rdurfey/R_misc/BayesianMethods/DBDA2Eprograms/DBDA2E-utilities.R')
suppressWarnings(library(rstan))

##### 2.1 Estimation of binomial probability

# Specify model:
modelString = "
data {
int<lower=0> N ;
int y[N] ; //length-N vector of integers 
}
parameters {
real<lower=0,upper=1> theta ;
}
model {
theta ~ beta(1,1) ;
y ~ bernoulli(theta) ; 
}
"
# notice above we use 'beta()' instead of 'dbeta()'


# Translate model to C++ and compile to DSO:
stanDso <- stan_model( model_code=modelString ) 

# Data are specified similar to JAGS.
# Specify data:
N = 50 ; z = 10
y = c(rep(1,z),rep(0,N-z))
dataList = list(
  y = y ,
  N = N 
)

# Running MCMC is done by sampling().
# Argument warmup has the same meaning as “burn in”.
# Argument ’iter` is total number of steps per chain, including warmup period.
# Argument init is not used in the following call allowing Stan to used default random initiation.
# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 3 ,
                     iter = 1000 , 
                     warmup = 200 , 
                     thin = 1 )

class(stanFit)

# Use application shinystan for exploration of the MCMC object.
suppressWarnings(library(shinystan))

#library(brms) # not sure what this is...

# Create an object of shinystan application.
# Save it.
myFirstStanFit<-as.shinystan(stanFit)
save(myFirstStanFit,file=paste(datapath,"./firststanfit.Rdata",sep="/"))

# Launch the application.
launch_shinystan(myFirstStanFit)

# Or load the file “firststanfit.Rdata” and launch the application later.
load(paste(datapath,"firststanfit.Rdata",sep="/"))
launch_shinystan(myFirstStanFit)

# Or use standard graphs from Stan or from [K].
traceplot(stanFit,pars=c("theta"))
#saveGraph(file=paste0(fileNameRoot,"StanTrace"),type="eps")
#openGraph()
plot(stanFit,pars=c("theta"))
#saveGraph(file=paste0(fileNameRoot,"StanPlot"),type="eps")

# Make graphs:
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format. This excludes warmup and thinned steps.
mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , 
                              function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
diagMCMC( mcmcCoda , parName=c("theta") )
#saveGraph(file=paste0(fileNameRoot,"Diag"),type="eps")



##### 2.2. Repeated use of the same DSO

# Make another data set.
# Specify data:
N = 50 ; z = 40
y = c(rep(1,z),rep(0,N-z))
dataList = list(
  y = y ,
  N = N 
)
dataList

# Run MCMC with these data.
# Note that we use the same dynamic shared object (DSO)
# Generate posterior sample:
stanFit <- sampling( object=stanDso , 
                     data = dataList , 
                     chains = 3 ,
                     iter = 1000 , 
                     warmup = 200 , 
                     thin = 1 )


# Use shinystan.
mySecondStanFit<-as.shinystan(stanFit)
save(mySecondStanFit,file=paste(datapath,"secondstanfit.Rdata",sep="/"))
launch_shinystan(mySecondStanFit)
load(paste(datapath,"secondstanfit.Rdata",sep="/"))
launch_shinystan(mySecondStanFit)


# Explore the graphs.
#openGraph()
traceplot(stanFit,pars=c("theta"))
#openGraph()
plot(stanFit,pars=c("theta"))
# Make graphs:
# For consistency with JAGS-oriented functions in DBDA2E collection, 
# convert stan format to coda format. This excludes warmup and thinned steps.
mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , 
                              function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
diagMCMC( mcmcCoda , parName=c("theta") )



##### 2.3. General structure of model description in Stan
...

##### 2.4. The following example shows how to sample from prior distribution


##### 3. Therapeutic touch example [K], page 240

source(paste(datapath,"./DBDA2Eprograms/Stan-Ydich-XnomSsubj-MbernBetaOmegaKappa.R",sep="/"))
show(genMCMC)