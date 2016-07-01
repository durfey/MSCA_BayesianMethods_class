################################
# Week 7 Workshop Notes & Code #
################################



# 1. Fitting normal model for 1 group with no predictors
# This example is based on [K], Chapter 16.

# preparing data
set.seed(9384756)
y <- rnorm(100, mean = 6, sd = 3)
Ntotal = length(y)

dataList = list(
  y = y ,
  Ntotal = Ntotal ,
  meanY = mean(y) ,
  sdY = sd(y)
)

source("//Users/rdurfey/R_Misc/BayesianMethods/DBDA2Eprograms/DBDA2E-utilities.R")


# 1.1. Running in JAGS
library(rjags)

# 
# modelString = "
# model {
# for(i in 1:Ntotal){
# y[i] ~ dnorm(mu,1/sigma^2)
# }
# mu ~ dnorm(mean=M,sd=S)
# sigma ~ dgamma(shape=alpha, rate=beta)
# }
# "


modelString = "
model {
for ( i in 1:Ntotal ) {
y[i] ~ dnorm( mu , 1/sigma^2 )
}
mu ~ dnorm( meanY , 1/(100*sdY)^2 )
sigma ~ dunif( sdY/1000 , sdY*1000 )
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )


# Initialize chains.
mu = mean(y)
sigma = sd(y) 
initsList = list( mu = mu , sigma = sigma)


# Run the chains.
parameters = c( "mu" , "sigma")     # The parameters to be monitored
adaptSteps = 500               # Number of steps to "tune" the samplers
burnInSteps = 1000
numSavedSteps=50000
nChains = 4 
thinSteps = 1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Run it
# The saved MCMC chain:
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )

# Check how chains converged.
summary(codaSamples)

# The parameters are esimated close to what we simulated and very similar to what point estimation would give.
mean(y)
sd(y)

# The plot of the samples and the densities of the parameters.
plot(codaSamples)

# The plot of autocorrelations.
autocorr.plot(codaSamples,ask=F)

# Autocorrelation function plot shows that standard deviation effective size must be pretty small.
effectiveSize(codaSamples)

# Shrink factor shows that even with long memory for standard deviation distributions converged:
gelman.diag(codaSamples)
gelman.plot(codaSamples)

(HDIofChains<-lapply(codaSamples,function(z) cbind(Mu=HDIofMCMC(codaSamples[[1]][,1]),Sd=HDIofMCMC(codaSamples[[1]][,2]))))



# 1.2. Running in Stan
library(rstan)


# Use the description for Stan from file “ch16_1.stan” or write it yourself and compare with the code.
modelString = "
data {
int<lower=1> Ntotal;
real y[Ntotal];
real meanY;
real sdY;
}
transformed data {
real unifLo;
real unifHi;
real normalSigma;
unifLo <- sdY/100;
unifHi <- sdY*100;
normalSigma <- sdY*100;
}
parameters {
real mu;
real<lower=0> sigma;
}
model {
sigma ~ uniform(unifLo, unifHi);
mu ~ normal(meanY, normalSigma); 
y ~ normal(mu, sigma);
}
" # close quote for modelString


stanDso <- stan_model( model_code=modelString ) 

stanFit <- sampling( object=stanDso , 
                     data = dataList ,
                     pars=c('mu', 'sigma'),
                     chains = 2,
                     cores= 2,
                     iter = 5000,
                     warmup = 200, 
                     thin = 1 )


# Alternatively, use the description in “ch16_1.stan”.
dataPath<-"//Users/rdurfey/R_misc/BayesianMethods/Week7"
# fit model
fit <- stan (file=paste(dataPath,"ch16_1.stan",sep="/"), 
             data=list(Ntotal=length(y),
                       y=y,
                       meanY=mean(y),
                       sdY=sd(y)), 
             pars=c('mu', 'sigma'),
             #control=list(adapt_delta=0.99),
             iter=5000, chains = 2, cores = 2
)

# Objects fit and stanFit should return very similar results.

# text statistics:
print (fit)

# estimates & hdi:
plot(fit)

# samples
traceplot(fit, ncol=1, inc_warmup=F)
pairs(fit, pars=c('mu','sigma'))
stan_scat(fit, c('mu', 'sigma'))
stan_hist(fit)
stan_dens(fit)

# autocorrelation:
stan_ac(fit, separate_chains = T)


# or you can work with familiar coda class:
library(coda)
stan2coda <- function(fit) {
  # apply to all chains
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
codaSamples <- stan2coda(fit)
summary(codaSamples)