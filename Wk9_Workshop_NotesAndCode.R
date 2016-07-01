################################
# Week 9 Workshop Notes & Code #
################################

##########################
# 1. Multiple Regression #
##########################

# Example 1

# Prepare data for multiple linear regression with 2 predictors.
# generate data
set.seed(3526)
Ntotal <- 1000
x <- cbind(rnorm(Ntotal, mean = 20, sd = 4), 
           rnorm(Ntotal, mean=10, sd = 6)
)
Nx <- ncol(x)
y <- 4 + .1*x[,1] + 3*x[,2] + rnorm(Ntotal, mean = 0, sd = 1)
dataListRegression<-list(Ntotal=Ntotal,
                         y=y,
                         x=as.matrix(x),
                         Nx=Nx)


library(rstan)

# Description of the model:
modelString<-"
data {
  int<lower=1> Ntotal;
  int<lower=1> Nx;
  vector[Ntotal] y;
  matrix[Ntotal, Nx] x;
}
transformed data {
  real meanY;
  real sdY;
  vector[Ntotal] zy; // normalized
  vector[Nx] meanX;
  vector[Nx] sdX;
  matrix[Ntotal, Nx] zx; // normalized
  
  meanY <- mean(y);
  sdY <- sd(y);
  zy <- (y - meanY) / sdY;
  for ( j in 1:Nx ) {
    meanX[j] <- mean(x[,j]);
    sdX[j] <- sd(x[,j]);
    for ( i in 1:Ntotal ) {
      zx[i,j] <- ( x[i,j] - meanX[j] ) / sdX[j];
    }
  }
}
parameters {
  real zbeta0;
  vector[Nx] zbeta;
  real<lower=0> nu;
  real<lower=0> zsigma;
}
transformed parameters{
  vector[Ntotal] zy_hat;
  zy_hat <- zbeta0 + zx * zbeta;
}
model {
  zbeta0 ~ normal(0, 2);
  zbeta  ~ normal(0, 2);
  nu ~ exponential(1/30.0);
  zsigma ~ uniform(1.0E-5 , 1.0E+1);
  zy ~ student_t(1+nu, zy_hat, zsigma);
}
generated quantities { 
  // Transform to original scale:
  real beta0; 
  vector[Nx] beta;
  real sigma;
  // .* and ./ are element-wise product and divide
  beta0 <- zbeta0*sdY  + meanY - sdY * sum( zbeta .* meanX ./ sdX );
  beta <- sdY * ( zbeta ./ sdX );
  sigma <- zsigma * sdY;
} "


stanDso<-stan_model( model_code=modelString )

# fit model
fit1<-sampling(stanDso,data=dataListRegression,pars=c('beta0', 'beta', 'nu', 'sigma'),iter=5000, chains = 2, cores = 2)

# Look at the results.
pairs(fit1)
plot(fit1)

# Analyze fitted model using shinystan
library(shinystan)
launch_shinystan(fit1)

summary(fit1)



# Example 2: insignificant predictor
# Read the data from homework assignment of week 7 of Statistical Analysis (31007).
# Run the chains and analyze the rezults.

Regression.Data<-as.matrix(read.csv("./Week9/DataForRegressionANOVA.csv"))
head(Regression.Data)

Ntotal <- nrow(Regression.Data)
x <- Regression.Data[,2:3]
head(x)

Nx <- ncol(x)
y <-Regression.Data[,1]
dataList<-list(Ntotal=Ntotal,
               y=y,
               x=as.matrix(x),
               Nx=Nx)

# Run MCMC.
fit2<-sampling(stanDso,data=dataList,pars=c('beta0', 'beta', 'nu', 'sigma'),iter=5000, chains = 2, cores = 2)
launch_shinystan(fit2)
pairs(fit2)
plot(fit2)
plot(fit2,pars=c('beta0', 'beta','sigma'))



# Example 3: Correlated predictors
# Create a data set with strongly correlated predictors.
set.seed(83945)
Ntotal <- 1000
x1 <- rnorm(Ntotal, mean = 20, sd = 4)
x2<-1-1.5*x1+rnorm(Ntotal, mean=0, sd = .1)
x<-cbind(x1,x2)           
plot(x)

Nx <- ncol(x)
y <- 4 + .2*x[,1] + 3*x[,2]+rnorm(Ntotal, mean = 0, sd = 1)
plot(x[,1],y)
plot(x[,2],y)

dataListShrink2<-list(Ntotal=Ntotal,
                      y=y,
                      x=as.matrix(x),
                      Nx=Nx)


# Run the chains and analyze the rezults.
tStart<-proc.time()
fit3<-sampling(stanDso,data=dataListShrink2,pars=c('beta0', 'beta', 'nu', 'sigma'),iter=5000, chains = 2, cores = 2)
tEnd<-proc.time()
tEnd-tStart

# Check convergence in shiny.
launch_shinystan(fit3)
pairs(fit3)
plot(fit3)
plot(fit3,pars=c('beta0', 'beta','sigma'))
stan_dens(fit3)
stan_ac(fit3, separate_chains = T)

betas<-cbind(Beta1=rstan::extract(fit3,pars="beta[1]")$'beta[1]',Beta2=rstan::extract(fit3,pars="beta[2]")$'beta[2]')

(lmFit<-lm(y~x1+x2))
summary(lmFit)
drop1(lmFit)




###########################################
# 2. Shrinkage of regression coefficients #
###########################################
# this is similar to Regularization techniques: Lasso/Ridge

# Use the same data dataListRegression as in the previous section.
modelString<-"
data {
int<lower=1> Ntotal;
int<lower=1> Nx;
vector[Ntotal] y;
matrix[Ntotal, Nx] x;
}
transformed data {
real meanY;
real sdY;
vector[Ntotal] zy; // normalized
vector[Nx] meanX;
vector[Nx] sdX;
matrix[Ntotal, Nx] zx; // normalized

meanY <- mean(y);
sdY <- sd(y);
zy <- (y - meanY) / sdY;
for ( j in 1:Nx ) {
meanX[j] <- mean(x[,j]);
sdX[j] <- sd(x[,j]);
for ( i in 1:Ntotal ) {
zx[i,j] <- ( x[i,j] - meanX[j] ) / sdX[j];
}
}
}
parameters {
real zbeta0;
real<lower=0> sigmaBeta;
vector[Nx] zbeta;
real<lower=0> nu;
real<lower=0> zsigma;
}
transformed parameters{
vector[Ntotal] zy_hat;
zy_hat <- zbeta0 + zx * zbeta;
}
model {
zbeta0 ~ normal(0, 2);
sigmaBeta ~ gamma(1,2); // mode 1.0, sd 1.0
zbeta  ~ student_t(1, 0, sigmaBeta);
nu ~ exponential(1/30.0);
zsigma ~ uniform(1.0E-5 , 1.0E+1);
zy ~ student_t(1+nu, zy_hat, zsigma);
}
generated quantities { 
// Transform to original scale:
real beta0; 
vector[Nx] beta;
real sigma;
// .* and ./ are element-wise product and divide
beta0 <- zbeta0*sdY  + meanY - sdY * sum( zbeta .* meanX ./ sdX );
beta <- sdY * ( zbeta ./ sdX );
sigma <- zsigma * sdY;
} "

stanShrinkDso<-stan_model( model_code=modelString )

tStart<-proc.time()
# fit model
fit4 <- sampling (stanShrinkDso, 
                  data=dataListShrink2, 
                  pars=c('beta0', 'beta', 'nu', 'sigma', 'sigmaBeta'),
                  iter=5000, chains = 2, cores = 2
)
tEnd<-proc.time()
tEnd-tStart


# Analyze fitted model using shinystan
launch_shinystan(fit4)
pairs(fit4)
plot(fit4,pars=c('beta0', 'beta','sigma','sigmaBeta'))
plot(fit3,pars=c('beta0', 'beta','sigma'))
stan_dens(fit3)
stan_ac(fit3, separate_chains = T)

summary(fit4)$summary
summary(fit3)$summary


