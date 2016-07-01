#################################
# Week 11 Workshop Notes & Code #
#################################

library(rstan)
source('//Users/rdurfey/R_misc/BayesianMethods/DBDA2Eprograms/DBDA2E-utilities.R')  

##############################
# 1. ANOVA in Bayesian Setup #
##############################

modelString<-"
data {
int<lower=1> Ntotal;
vector[Ntotal] y;
int<lower=2> Nx1Lvl;
int<lower=2> Nx2Lvl;
int<lower=1, upper=Nx1Lvl> x1[Ntotal];
int<lower=1, upper=Nx2Lvl> x2[Ntotal];
real<lower=0> agammaShRa[2];
}
transformed data {
real meanY;
real sdY;
vector[Ntotal] zy;
meanY <- mean(y);
sdY <- sd(y);
zy <- (y - mean(y)) / sdY;  // center & normalize
}
parameters {
real a0;
real<lower=0> a1Sigma;
real<lower=0> a2Sigma;
real<lower=0> a1a2Sigma;
vector[Nx1Lvl] a1;
vector[Nx2Lvl] a2;
matrix[Nx1Lvl,Nx2Lvl] a1a2;
real<lower=0> zySigma;
}
model {
a0 ~ normal(0, 1);
a1Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
a1 ~ normal(0, a1Sigma);
a2Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
a2 ~ normal(0, a2Sigma);
a1a2Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
for (j1 in 1:Nx1Lvl) {
a1a2[j1,] ~ normal(0, a1a2Sigma);
}
zySigma ~ uniform(1.0/10, 10);
for ( i in 1:Ntotal ) {
zy[i] ~ normal(a0 + a1[x1[i]] + a2[x2[i]]+ a1a2[x1[i],x2[i]], zySigma);
}
}
generated quantities {
// Convert a to sum-to-zero b :
real b0;
vector[Nx1Lvl] b1;
vector[Nx2Lvl] b2;
matrix[Nx1Lvl,Nx2Lvl] b1b2;
matrix[Nx1Lvl,Nx2Lvl] m;
real<lower=0> b1Sigma;
real<lower=0> b2Sigma;
real<lower=0> b1b2Sigma;
real<lower=0> ySigma;
for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2]; // cell means 
} }
b0 <- mean(m);
for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( m[j1,] ) - b0; }
for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( m[,j2] ) - b0; }
for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] );  
} }
// transform to original scale:
b0 <- meanY + sdY * b0;
b1 <- sdY * b1;
b2 <- sdY * b2;
b1b2 <- sdY * b1b2;
b1Sigma <- sdY * a1Sigma;
b2Sigma <- sdY * a2Sigma;
b1b2Sigma <- sdY * a1a2Sigma;
ySigma <- sdY * zySigma;
}"

# create DSO
stanDsoANOVA2Way<-stan_model( model_code=modelString )



# 20_1: Metric Predicted Variable with Two Nominal Predictors

# load data from 'Salary.csv' (see Kruschke)
mydf = read.csv("//Users/rdurfey/R_misc/BayesianMethods/Week11/Salary.csv")
mean(mydf$Salary)
head(mydf)
dim(mydf)
colnames(mydf)
table(mydf$Pos)
table(mydf$Org)
length(table(mydf$Org))

# the output will be salary
y <- mydf$Salary;
x1 <- mydf$Pos;
x2 <- mydf$Org;
dataListSalary<-list(Ntotal=length(y),
                     y=y,
                     x1=as.integer(x1),
                     x2=as.integer(x2),
                     Nx1Lvl=nlevels(x1),
                     Nx2Lvl=nlevels(x2),
                     agammaShRa=unlist( gammaShRaFromModeSD(mode=1/2, sd=2) ))


# Create names of variables and their interactions for further reference.
namesPos<-names(table(mydf$Pos))
namesOrg<-names(table(mydf$Org))
as.vector(outer(1:4,1:2,paste,sep="-"))

namesInter<-as.vector(outer(namesOrg,namesPos,paste,sep="-"))
varNames<-c("Intercept",namesPos,namesOrg,namesInter,rep("Var",5))

# Run MCMC
# fit model
fit <- sampling (stanDsoANOVA2Way, 
                 data=dataListSalary, 
                 pars=c('b0', 
                        'b1', 
                        'b2', 
                        'b1b2',
                        'b1Sigma', 
                        'b2Sigma',
                        'b1b2Sigma',
                        'ySigma'),
                 iter=5000, chains = 2, cores = 2
)

# Check the results in shinystan.
library(shinystan)
launch_shinystan(fit)

# Create results including mean value, 2.5%, 50% and 97.5% quantiles.
# Add variable names as row names.
SalaryResults<-summary(fit)$summary[,c(1,4,6,8)]
varNames[nrow(SalaryResults)-(4:0)]<-rownames(SalaryResults)[nrow(SalaryResults)-(4:0)]
rownames(SalaryResults)<-varNames
SalaryResults

plot(fit,pars=c("b1"))
plot(fit,pars=c('b2'))
plot(fit,pars=c("b1b2"))


# Extract chains for the position variables.
fit_ext <- rstan::extract(fit)
names(fit_ext)

fit_ext.b1<-fit_ext$b1
colnames(fit_ext.b1)<-namesPos
head(fit_ext.b1)


# Extract chains for the department variables.
fit_ext.b2<-fit_ext$b2
colnames(fit_ext.b2)<-namesOrg
head(fit_ext.b2)


# Extract chains for interaction variables.
fit_ext.b1.b2<-fit_ext$b1b2
dim(fit_ext.b1.b2)

dimnames(fit_ext.b1.b2)[[2]]<-namesPos
dimnames(fit_ext.b1.b2)[[3]]<-namesOrg
dimnames(fit_ext.b1.b2)

fit_ext.b1.b2[1,,]


#############
# EXERCISES #
#############

# 1. Use contrasts to compare salaries at Business and Finance with Physics and with Chemistry departments.
contrast_BFIN_PHYS <- fit_ext.b1.b2[,,"BFIN"] - fit_ext.b1.b2[,,"PHYS"]
plot(contrast_BFIN_PHYS)
hist(contrast_BFIN_PHYS)

contrast_BFIN_CHEM <- fit_ext.b1.b2[,,"BFIN"] - fit_ext.b1.b2[,,"CHEM"]
plot(contrast_BFIN_CHEM)
hist(contrast_BFIN_CHEM)


# 2. Use contrasts to compare salaries of Endowment full Professor and Distinguished Full professor.
contrast_NDW_DST <- fit_ext.b1.b2[,"NDW",] - fit_ext.b1.b2[,"DST",]
plot(contrast_NDW_DST)
hist(contrast_NDW_DST)
# conclusion: they make about the same salaries


# 3. Use contrasts to compare salaries spreads between Full Professor and Assistant Professor at Physics Department and at Chemistry Department.
contrast_FT1_FT3_PHYS <- fit_ext.b1.b2[,"FT1","PHYS"] - fit_ext.b1.b2[,"FT3","PHYS"]
plot(contrast_FT1_FT3_PHYS)
hist(contrast_FT1_FT3_PHYS)

contrast_FT1_FT3_CHEM <- fit_ext.b1.b2[,"FT1","CHEM"] - fit_ext.b1.b2[,"FT3","CHEM"]
plot(contrast_FT1_FT3_CHEM)
hist(contrast_FT1_FT3_CHEM)
# note: we think FT1 is full prof and FT3 is assist. prof
# conclusion: full prof makes more. duh.

# 4. Analyze contrasts for comparison of salary spreads between the departments of Physics and Chemistry.
contrast_contrasts<-contrast_FT1_FT3_PHYS - contrast_FT1_FT3_CHEM
hist(contrast_contrasts)
# conclusion: CHEM has a higher spread of salaries

##############################################################################
# 2. Understanding the effect of scaling and transformations on interactions #
##############################################################################


# Nonlinear transformations may affect interactions very significantly.

# Illustrate it on a simple simulated example.
mean00<-1
mean10<-3
mean01<-4
mean11<-6
y00<-rnorm(5,mean00,.1)
y10<-rnorm(5,mean10,.1)
y01<-rnorm(5,mean01,.1)
y11<-rnorm(5,mean11,.1)


# Plot the effects. If the lines are parallel the effects are additive.
plot(c(0,1),c(mean(y00),mean(y10)),type="b",ylim=c(1,8),col="darkgreen",lwd=3,ylab="Response",xlab="Predictor 1")
lines(c(0,1),c(mean(y01),mean(y11)),type="b",col="lightblue",lwd=3)
legend("topleft",legend=c("Predictor2 at 0","Predictor2 at 1"),lty=1,lwd=3,col=c("darkgreen","lightblue"))

# plot shows no interaction



# Taking exponent of the same data introduces significant interaction.
plot(c(0,1),c(mean(exp(y00)),mean(exp(y10))),type="b",ylim=c(1,400),col="darkgreen",lwd=3,ylab="Response",xlab="Predictor 1")
lines(c(0,1),c(mean(exp(y01)),mean(exp(y11))),type="b",col="lightblue",lwd=3)
legend("topleft",legend=c("Predictor2 at 0","Predictor2 at 1"),lty=1,lwd=3,col=c("darkgreen","lightblue"))

# now, plot shows significant interaction. the slope of the response depends on the level of the 2nd predictor

