##################################
# Week 2 Workshop - Notes & Code #
##################################

#-----------------------------------------------------
# MAMMOGRAPHY SCREENING CONTROVERSRY
# According to CDC, probability of breast cancer among women of age 50-60 years old is 2.28%.
# Screening using mammography has sensitivity (i.e. probability of detecting the disease correctly) 75%-94%. Specificity of the test (i.e. probability of correct “no disease” conclusion) is 83%-98%. Assume the most optimistic parameters of the test (94% and 98%, correspondingly for sensitivity and specificity) What is the probability that randomly selected woman with positive test result has the disease?
#-----------------------------------------------------


#Prior p(theta1)
p.theta1<-.0228
#Sensitivity: p(tau1|theta1)
p.tau1.theta1<-.94
#Specificity: p(tau0|theta0)
p.tau0.theta0<-.98
#Complimentary probabilities
p.theta0<-1-p.theta1
p.tau1.theta0<-1-p.tau0.theta0
#Bayes
(p.theta1.tau1<-p.tau1.theta1*p.theta1/(p.tau1.theta1*p.theta1+p.tau1.theta0*p.theta0))


#-------------------------------------------------
# Assume that randomly selected woman tested positive and then retested with 
# negative test result.
# After both tests what is the probability of disease?
#-------------------------------------------------


# Set the prior to the new probability of having the disease:
pDisease = p.theta1.tau1
# Bayes rule for second, negative test:
(p.theta1.tau0 = ((1.0-p.tau1.theta1) * pDisease / 
                    ( (1.0-p.tau1.theta1) * pDisease + (1.0-p.tau1.theta0) *
                          (1.0-pDisease))))
# NOTE: (1.0-p.tau1.theta1) is the equivalent of just p.tau0.theta1. 
# They are complements. Just didn't create a new variable above.


# MY SIDENOTE: what if she tested Negative, and then Positive?
p.tau0.theta1<-1-p.tau1.theta1
(pDisease2 <- p.tau0.theta1*p.theta1 / (p.tau0.theta0*p.theta0 + p.tau0.theta1*p.theta1))
(p.theta1.tau1_2 <- ((p.tau1.theta1)*pDisease2 / (p.tau1.theta1*pDisease2 + p.tau1.theta0*(1-pDisease2))))
# It's the same (which makes sense)


#-------------------------------------
# COMPARISON WITH FNP STATISTICAL APPROACH
# Example from Linear and Non-Linear Models, Lecture 4
#--------------------------------------

library(faraway)
data(babyfood)
babyfood
(data<-babyfood[c(1,3),c(1:2,4)])

#------------------
# FNP Approach
#-----------------

prop.test(as.matrix(data[,1:2]))
# The null hypothesis is rejected.

# Then we use logistic regression to find out how much the odds of having the disease change between bottle-fed and breast-fed babies.
mdl<-glm(cbind(disease,nondisease)~food,family=binomial,data)
summary(mdl)

# Take exponents to express the change in odds ratio, rather than log of it.
exp(mdl$coefficients)
# NOTE: the second coefficient is the multiplicative change from the intercept
# Here, this is saying that the odds of a breast-fed baby boy having the disease is 
# approx. HALF of the odds shown by the intercept (which is bottle fed

# Look at confidence intervals to see the significance and the range of accuracy of our conclusions.
library(MASS)
exp(confint(mdl))
# NOTE: it's important to note that the foodBreast interval does NOT include 1 (one),
# because that would indicate that there is a chance that breast-feeding could INCREASE
# odds of disease. That would basically invalidate the conclusion and we couldn't reject
# the null hypothesis that there isn't a difference between food & breast.


#-------------------------
# Bayesian Approach
# BUT, it's a not-quite-true bayesian approach..
#-------------------------

# Create joint distribution of two binary variables: (“disease”,“nondisease”) and (“bottle”,“breast”).
joint.dist<-data
joint.dist[1:2,1:2]<-joint.dist[1:2,1:2]/sum(joint.dist[1:2,1:2])
joint.dist

# Find marginal probabilities the disease and the treatment.
(p.breast<-sum(joint.dist[2,-3]))
(p.disease<-sum(joint.dist[,1]))

# Find conditional probability of breastfeeding, given that baby got the disease
(p.breast.disease<-joint.dist[2,1]/p.disease)

# Finally, use Bayes theorem.
(p.disiease.breast<-p.breast.disease*p.disease/p.breast)

# Could we calculate the same probability directly, not using Bayes theorem?
(p.disiease.breast<-joint.dist[2,1]/p.breast)

# QUESTION: Why one would prefer using Bayes theorem rather than calculating the probability directly?
# ANSWER: The difference (and preference) lies in how we collect the data..

# QUESTION: This is not Bayesian approach. Why?
# ANSWER: Because we're not calculating a prior distribution. We're just solving the contingency table using Bayes Theorum

#----------------------------------------------
# True B-approach but for a different problem
#----------------------------------------------

# The data collected from observation of 952 baby boys of which 124 got the disease: this is a sequence of 124 ones and 828 zeros.
(Data<-c(rep(0,828),rep(1,124)))

# Probability of success is p.
# In the FNP-approach its estimate is p̂.
# In the B-approach it is a random variable with prior distribution density fp(x)fp(x).
# Then the B-approach means using Bayes theorem to find posterior fp(p|Data)fp(p|Data).
# We will practice this approach in the following section.

#--------------------------------------------------
# Example of Using Bayes Theorem in Statistics
#--------------------------------------------------
dataPath<-"./"
source(paste(dataPath,"DBDA2Eprograms/DBDA2E-utilities.R",sep="/"))
source(paste(dataPath,"DBDA2Eprograms/BernGrid.R",sep="/"))

# Binomial model, triangular prior with 3 values

(Theta = seq( 0 , 1 , length=5 ))  # Sparse teeth for Theta.

pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
(pTheta = pTheta/sum(pTheta)) # this is the prior

Data = c(rep(0,0),rep(1,1)) 

# likelihood & posterior
(likelihood <- Theta^sum(Data)*(1-Theta)^(1-sum(Data)))
post<-pTheta*likelihood
(post<-post/sum(post))

# Use the BernGrid function (sourced above)
#openGraph(width=5,height=7)
(posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                       showCentTend="None" , showHDI=FALSE , showpD=FALSE ))


# Binomial model, TRIANGULAR prior with 11 values

Theta = seq( 0 , 1 , length=11 )  # Sparse teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))      # Single flip with 1 head

#openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )

# QUESTION: Which of the functions has more influence on the posterior: prior or likelihood?
# ANSWER: The posterior is still influenced more by the prior, but not entirely so.

# NOTE: notice that the Posterior is more affected by the Likelihood here than
# it was in the previous example.
# The larger the sample, the more the Likelihood will influence the posterior


# Binomial model, UNIFORM prior with 1001 values

Theta = seq( 0 , 1 , length=1001 ) # Fine teeth for Theta.
pTheta = rep(1,length(Theta))      # Uniform (horizontal) shape for pTheta.
pTheta = pTheta/sum(pTheta)        # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))        # Single flip with 1 head

posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )

# QUESTION: Which of the functions has more influence on the posterior: prior or likelihood?
# ANSWER: The Likelihood (obviously).


# Binomial model, binary prior with values 0, 1

Theta = seq( 0 , 1 , length=1001 ) # Fine teeth for Theta.
pTheta = rep(0,length(Theta))      # Only extremes are possible!
pTheta[2] = 1                      # Only extremes are possible!
pTheta[length(pTheta)-1] = 1       
pTheta = pTheta/sum(pTheta)        # Make pTheta sum to 1.0
Data = c(rep(0,0),rep(1,1))        # Single flip with 1 head

#openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=FALSE , showpD=FALSE )


# Binomial model, triangular prior with 1001 values, 4 observations

Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

# first, do it yourself... likelihood & posterior
likelihood <- Theta^sum(Data)*(1-Theta)^(4-sum(Data))
post<-pTheta*likelihood
post<-post/sum(post)
plot(Theta,pTheta)
plot(Theta,likelihood)
plot(Theta,post)

#openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="None" , showHDI=TRUE , showpD=FALSE )



# Binomial model, concentrated prior, 4 observations

Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^10               # Sharpen pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

#openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )




# Binomial model, dispersed prior, 4 observations

Theta = seq( 0 , 1 , length=1001 )  # Fine teeth for Theta.
pTheta = pmin( Theta , 1-Theta ) # Triangular shape for pTheta.
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
pTheta = pTheta^0.1              # Flatten pTheta !
pTheta = pTheta/sum(pTheta)      # Make pTheta sum to 1.0
Data = c(rep(0,3),rep(1,1))      # 25% heads, N=4

#openGraph(width=5,height=7)
posterior = BernGrid( Theta, pTheta , Data , plotType="Bars" , 
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )







#--------------------------------------------------------------------------
# Posterior distribution by simulation. Updating prior and adding new data
# Look again at the example in section 3.3. with uniform prior.
# Follow the steps of obtaining posterior distribution by Monte Carlo.
#--------------------------------------------------------------------------

# Define likelihood function for binomial distribution.
likeli<-function(par,data){
  sdata<-sum(data)
  ldata<-length(data)
  return(par^sdata*(1-par)^(ldata-sdata))
}

# Define values of parameter θ and prior distribution.
Theta = seq( .00001 , 1-.00001 , length=1001 ) # Fine teeth for Theta.
pTheta = rep(1,length(Theta))      # Uniform (horizontal) shape for pTheta.
pTheta = pTheta/sum(pTheta)        # Make pTheta sum to 1.0
plot(Theta,pTheta)

# Create data of length 5 generated by the model with parameter 0.84.
set.seed(5)
(Data<-rbinom(5,size=1,prob=.84))

# Create sample of θ generated from the prior distribution.
set.seed(15)
priorInd<-sample(1:length(Theta),500,replace = T)
priorSample<-cbind(Theta=Theta[priorInd],Prob=pTheta[priorInd])
priorSample<-rbind(priorSample,
                   c(head(Theta,1),head(pTheta,1)),
                   c(tail(Theta,1),tail(pTheta,1)))

# Calculate likelihood for each simulated θθ and the data.
likelihoodVector<-sapply(priorSample[,"Theta"],function(z) likeli(z,Data))
plot(priorSample[,"Theta"],likelihoodVector)

# Calculate posterior distribution.

# Calculate vector of numerators of the Bayes theorem
# Normalize it
# Create function for linear interpolation of vector of numerator
postVector<-priorSample[,"Prob"]*likelihoodVector
postVector<-postVector/sum(postVector)
plot(priorSample[,"Theta"],postVector)

postDistr<-approxfun(priorSample[,"Theta"],postVector,method="linear")
plot(priorSample[,"Theta"],postVector)
lines(Theta,postDistr(Theta),col="red",lwd=2)
