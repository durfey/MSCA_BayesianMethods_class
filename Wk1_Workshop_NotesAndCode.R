##################################
# Week 1 Workshop - Notes & Code #
##################################
dataPath<-"/Users/rdurfey/R_misc/BayesianMethods/Week1"

# Create a template for 3-dimensional distribution table.
u<-c("u1","u2")
v<-c("v1","v2","v3")
w<-c("w1","w2","w3")
matr.u0<-paste("u0",outer(v,w,paste,sep=","),sep=",")
dim(matr.u0)<-c(3,3)
matr.u1<-paste("u1",outer(v,w,paste,sep=","),sep=",")
dim(matr.u1)<-c(3,3)
matr.u0
matr.u1

data3way<-read.csv(file=paste(dataPath,"3WayData.csv",sep="/"))
head(data3way)

mat.u0<-table(subset(data3way,u==0)[,1],subset(data3way,u==0)[,2])
mat.u1<-table(subset(data3way,u==1)[,1],subset(data3way,u==1)[,2])
mat.u0
mat.u1

# Check elements of mat.u1 and name the columns and rows of both matrices.
idx.v1<-data3way$v==1
idx.w1<-data3way$w==1
idx.u1<-data3way$u==1
sum(idx.v1*idx.w1*idx.u1) #element (1,1) of mat.u1

idx.v2<-data3way$v==2
sum(idx.v2*idx.w1*idx.u1) #element (1,2) of mat.u1

idx.w2<-data3way$w==2
sum(idx.v1*idx.w2*idx.u1) #element (2,1) of mat.u1

colnames(mat.u1)<-colnames(mat.u0)<-c("v1","v2","v3")
rownames(mat.u1)<-rownames(mat.u0)<-c("w1","w2","w3")

(data3way.table<-list(u0=mat.u0,u1=mat.u1))

# Create 3-dimensional joint distribution.
N<-sum(unlist(data3way.table))
(data3way.table.p<-lapply(data3way.table,"/",N)) # counts divided by the total


#################################
# Create Marginal Distributions

# u
uMarginal<-c(u0=sum(data3way.table.p$u0),u1=sum(data3way.table.p$u1))
uMarginal

# v
marg.v1<-sum(data3way.table.p$u0[,1],data3way.table.p$u1[,1])
marg.v2<-sum(data3way.table.p$u0[,2],data3way.table.p$u1[,2])
marg.v3<-sum(data3way.table.p$u0[,3],data3way.table.p$u1[,3])
vMarginal<-c(v1=marg.v1,v2=marg.v2,v3=marg.v3)
vMarginal

# w
marg.w1<-sum(data3way.table.p$u0[1,],data3way.table.p$u1[1,])
marg.w2<-sum(data3way.table.p$u0[2,],data3way.table.p$u1[2,])
marg.w3<-sum(data3way.table.p$u0[3,],data3way.table.p$u1[3,])
wMarginal<-c(w1=marg.w1,w2=marg.w2,w3=marg.w3)
wMarginal

############################
# Create Conditional Distributions

# Create conditional distribution  p(w,v|u=1)
(cond.v.w.given.u1<-data3way.table.p[["u1"]])

# Create conditional distribution p(v|u=1)
(cond.v.given.u1<-apply(data3way.table.p[["u1"]],2,sum)/uMarginal["u1"])

sum(cond.v.given.u1)

# Create conditional distribution p(w|v=2,u=1) ----> p(w,v,u) = p(w|v,u)p(v,u) --> p(w|v,u)p(v|u)p(u)
(cond.w.given.u1.v2<-(data3way.table.p[["u1"]][,"v2"])/(cond.v.given.u1[["v2"]]*uMarginal["u1"]))

# Compare the vectors p(w|v2,u1)p(v2|u1)p(u1)p(w|v2,u1)p(v2|u1)p(u1) and p(w,v,u)[,v2,u1]
rbind(uMarginal["u1"]*cond.v.given.u1["v2"]*cond.w.given.u1.v2,data3way.table.p[["u1"]][,"v2"])

#-------------------------------------------------#
# MORAL OF THE STORY: p(w,v,u)=p(w|v,u)p(v|u)p(u) #
#-------------------------------------------------#

# similarly...
# p(u|v,w)=p(u,v,w)/p(v|w)p(w)


##############################################
# Simulation Using Conditional Distributions #
##############################################

# Let the marginal distribution for random variable vv be Bernoulli with p(u=0)=0.55p(u=0)=0.55.
# Let conditional distributions for random variables (v|u=0)(v|u=0) and (v|u=1)(v|u=1), taking values 1,2,31,2,3 be

(pCond.v.given.u0<-c(.7,.2,.1))
(pCond.v.given.u1<-c(.1,.2,.7))

# Let random variable (w|v,u)(w|v,u) take values 1,2,31,2,3 with probabilities p(w|v,u)p(w|v,u), given by the following:
p.given.u0.v1<-c(.3,.3,.4)
p.given.u0.v2<-c(.5,.3,.2)
p.given.u0.v3<-c(.6,.2,.2)
p.given.u1.v1<-c(.2,.3,.5)
p.given.u1.v2<-c(.2,.2,.6)
p.given.u1.v3<-c(.1,.7,.2)

# Simulate joint sample (w,v,u)(w,v,u) of lenth n=500n=500.
# Use set.seed(11) Start with simulation of uu.
# For each simulated value uu generate vv from the corresponding distribution.
# Finally, for each pair v,uv,u simulate ww.

set.seed(11)

sim.u<-rbinom(500,1,0.55)
sim.v<-rbinom(500,2,pCond.v.given.u0+pCond.v.given.u1)
sim.w<-
  
  pCond.v.given.u0+pCond.v.given.u1
