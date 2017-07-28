#The model code below is written for program R and uses the R2WinBUGS package 
#to run WinBUGS as well as the reshape package to format the occurrence data.

#It is designed to estimate static species-specific occupancy and detection 
#with site specific habitat and sampling covariates using the community model.
#The occurence data are in the file "occ data.csv". 
#The covariate data are in the files "habitat.csv" (occurence) and "date.csv" (detection).
#Species are grouped into one of three categories in the file "groups.csv"
#See Zipkin et al. 2010 (Biological Conservation) for more context and details on the model. 

#Read in the occurence data
data1 <- read.table("occ data.csv", header=TRUE,sep=",",na.strings=TRUE)
data1$Occ <- rep(1, dim(data1)[1])
#See the first ten lines of data
data1[1:10,]
#How many citings for each species
total.count = tapply(data1$Occ, data1$Species, sum)

#Find the number of unique species
uspecies = as.character(unique(data1$Species))
#n is the number of observed species
n=length(uspecies)

#Load the species groups data
groups <- read.table("groups.csv", header=TRUE,sep=",",na.strings=c("NA"))
species = as.character(groups$species)
assmb = groups$group
a=which(assmb==1)
ground = assmb; ground[-a] = 0
b=which(assmb==2)
mid=assmb; mid[-b] = 0; mid[b]=1

#Find the number of unique sampling locations
upoints = as.character(unique(data1$Point))
#J is the number of sampled points
J=length(upoints)

#Reshape the data using the R package "reshape"
library(reshape)

#The detection/non-detection data is reshaped into a three dimensional 
#array X where the first dimension, j, is the point; the second 
#dimension, k, is the rep; and the last dimension, i, is the species. 
junk.melt=melt(data1,id.var=c("Species", "Point", "Rep"), measure.var="Occ")
X=cast(junk.melt, Point ~ Rep ~ Species)

#Add in the missing lines with NAs
for (i in 1: dim(X)[3]) {
   b = which(X[,,i] > 0) 
   X[,,i][b] = 1  
   X[,,i][-b] = 0  
   X[,,i][1:36,4] = NA;  X[,,i][38:56,4] = NA;  
   X[,,i][59:61,4] = NA;  X[,,i][66:70,4] = NA;        
}

#Create all zero encounter histories to add to the detection array X 
#as part of the data augmentation to account for additional 
#species (beyond the n observed species). 

#nzeroes is the number of all zero encounter histories to be added
  nzeroes = 50
#X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled  
  X.zero = matrix(0, nrow=70, ncol=4)
  X.zero[1:36,4] = NA;  X.zero[38:56,4] = NA;  
  X.zero[59:61,4] = NA;  X.zero[66:70,4] = NA;   
#Xaug is the augmented version of X.  The first n species were actually observed
#and the n+1 through nzeroes species are all zero encounter histories  
  Xaug <- array(0, dim=c(dim(X)[1],dim(X)[2],dim(X)[3]+nzeroes))
  Xaug[,,(dim(X)[3]+1):dim(Xaug)[3]] = rep(X.zero, nzeroes)
  dimnames(X)=NULL
  Xaug[,,1:dim(X)[3]] <-  X

#K is a vector of length J indicating the number of reps at each point j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)

#Create a vector to indicate which habitat type each point is in (CATO = 1; FCW =0)
Ind <- as.vector(cbind(matrix(rep(0,35),ncol=1,nrow=35),
      matrix(rep(1,35),ncol=1,nrow=35)));

#Read in the habitat data      
habitat <- read.table("habitat.csv", header=TRUE,sep=",",na.strings=c("NA"))

#Standardize the understory foliage data (ufc)
ufc <- as.vector(habitat$ufc)
mufc <- mean(ufc, na.rm=TRUE)
sdufc <- sd(ufc, na.rm=TRUE)
ufc1 <- as.vector( (ufc-mufc) / sdufc )
ufc2 <- as.vector( ufc1*ufc1 )

#Standardize the tree basal area data (ba)
ba <- as.vector(habitat$ba)
mba <- mean(ba, na.rm=TRUE)
sdba <- sd(ba, na.rm=TRUE)
ba1 <- as.vector( (ba-mba) / sdba )
ba2 <- as.vector( ba1*ba1 )

#Read in the date data
#The sampling dates have been converted to Julien dates
dates <- read.table("dates.csv", header=TRUE,sep=",",na.strings=c("NA"))
dates <- as.matrix(dates[,c("rep1","rep2","rep3","rep3")])
mdate <- mean(dates, na.rm=TRUE)
sddate <- sqrt(var(dates[1:length(dates)], na.rm=TRUE))
date1 <- (dates-mdate) /  sddate
date2 <- date1*date1
date1 <- as.matrix(date1)
date2 <- as.matrix(date2)

#Write the model code to a text file 
cat("
   model{

#Define prior distributions for community-level model parameters
omega ~ dunif(0,1)

cato.mean ~ dunif(0,1)
mu.ucato <- log(cato.mean) - log(1-cato.mean)

fcw.mean ~ dunif(0,1)
mu.ufcw <- log(fcw.mean) - log(1-fcw.mean)

cato2.mean ~ dunif(0,1)
mu.vcato <- log(cato2.mean) - log(1-cato2.mean)

fcw2.mean ~ dunif(0,1)
mu.vfcw <- log(fcw2.mean) - log(1-fcw2.mean)

mua1 ~ dnorm(0, 0.001)
mua2 ~ dnorm(0, 0.001)
mua3 ~ dnorm(0, 0.001)
mua4 ~ dnorm(0, 0.001)
mub1 ~ dnorm(0, 0.001)
mub2 ~ dnorm(0, 0.001)

tau.ucato ~ dgamma(0.1,0.1)  
tau.ufcw ~ dgamma(0.1,0.1)
tau.vcato ~ dgamma(0.1,0.1) 
tau.vfcw ~ dgamma(0.1,0.1)
tau.a1 ~ dgamma(0.1,0.1)
tau.a2 ~ dgamma(0.1,0.1)
tau.a3 ~ dgamma(0.1,0.1)
tau.a4 ~ dgamma(0.1,0.1) 
tau.b1 ~ dgamma(0.1,0.1) 
tau.b2 ~ dgamma(0.1,0.1)

for (i in 1:(n+nzeroes)) {

#Create priors for species i from the community level prior distributions
    w[i] ~ dbern(omega)
    u.cato[i] ~ dnorm(mu.ucato, tau.ucato)
    u.fcw[i] ~ dnorm(mu.ufcw, tau.ufcw)  
    v.cato[i] ~ dnorm(mu.vcato, tau.vcato) 
    v.fcw[i] ~ dnorm(mu.vfcw, tau.vfcw)   
    a1[i] ~ dnorm(mua1, tau.a1)
    a2[i] ~ dnorm(mua2, tau.a2)
    a3[i] ~ dnorm(mua3, tau.a3)
    a4[i] ~ dnorm(mua4, tau.a4)     
    b1[i] ~ dnorm(mub1, tau.b1)    
    b2[i] ~ dnorm(mub2, tau.b2)


#Create a loop to estimate the Z matrix (true occurrence for species i 
#at point j.      
   for (j in 1:J) {
       logit(psi[j,i]) <- u.cato[i]*(1-Ind[j]) + u.fcw[i]*Ind[j] + 
               a1[i]*ufc1[j] + a2[i]*ufc2[j] + a3[i]*ba1[j] + a4[i]*ba2[j] 
       
  mu.psi[j,i] <- psi[j,i]*w[i]
  Z[j,i] ~ dbern(mu.psi[j,i])

#Create a loop to estimate detection for species i at point k during 
#sampling period k.      
     for (k in 1:K[j]) {  
    logit(p[j,k,i]) <-  v.cato[i]*(1-Ind[j]) + v.fcw[i]*Ind[j] + 
                      b1[i]*date1[j,k] + b2[i]*date2[j,k] 

       mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
       X[j,k,i] ~ dbern(mu.p[j,k,i])
}   }}


#Sum all species observed (n) and unobserved species (n0) to find the 
#total estimated richness
n0 <- sum(w[(n+1):(n+nzeroes)])
N <- n + n0


#Create a loop to determine point level richness estimates for the 
#whole community and for subsets or assemblages of interest.
for(j in 1:J){
Nsite[j]<- inprod(Z[j,1:(n+nzeroes)],w[1:(n+nzeroes)])
Nground[j]<- inprod(Z[j,1:n],ground[1:n])
Nmid[j]<- inprod(Z[j,1:n],mid[1:n])
}

#Finish writing the text file into a document we call covarmodel.txt
}
",file="covarmodel.txt")


#Load the R2Winbugs library
library(R2WinBUGS)

#Create the necessary arguments to run the bugs() command 
#Load all the data
sp.data = list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug, date1=date1,
               date2=date2, ufc1=ufc1, ba1=ba1, ufc2=ufc2, 
               ba2=ba2, Ind=Ind, ground=ground, mid=mid)

#Specify the parameters to be monitored
sp.params = list('u.cato', 'u.fcw', 'v.cato', 'v.fcw', 'omega', 'a1', 
		'a2', 'a3', 'a4', 'b1', 'b2', 'Nsite', 'N', 'Nground', 'Nmid') 

#Specify the initial values
    sp.inits = function() {
    omegaGuess = runif(1, n/(n+nzeroes), 1)
    psi.meanGuess = runif(1, .25,1)
    list(omega=omegaGuess,w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
               u.cato=rnorm(n+nzeroes), v.cato=rnorm(n+nzeroes),
               u.fcw=rnorm(n+nzeroes), v.fcw=rnorm(n+nzeroes),
               Z = matrix(rbinom((n+nzeroes)*J, size=1, prob=psi.meanGuess), 
		                nrow=J, ncol=(n+nzeroes)), 
               a1=rnorm(n+nzeroes), a2=rnorm(n+nzeroes), a3=rnorm(n+nzeroes), 
               a4=rnorm(n+nzeroes), b1=rnorm(n+nzeroes), b2=rnorm(n+nzeroes)
               )
           }
 
#Run the model and call the results “fit”
fit = bugs(sp.data, sp.inits, sp.params, "covarmodel.txt", debug=TRUE, 
         n.chains=2, n.iter=10000, n.burnin=5000, n.thin=5)

#######################################################################
#Summarize some results

#See a summary of the parameter estimates
fit$summary

#See baseline estimates of species-specific occupancy and detection in one of 
#the habitat types (CATO)
cato.occ = fit$sims.list$u.cato
cato.det = fit$sims.list$v.cato

#This includes occupancy and detection estimates for all observed 
#species only (species 1:n)
psi.cato = plogis(cato.occ[,1:n]) 
p.cato   = plogis(cato.det[,1:n]) 

occ.matrix <- cbind(apply(psi.cato,2,mean),apply(psi.cato,2,sd))
det.matrix <- cbind(apply(p.cato,2,mean),apply(p.cato,2,sd))

#See estimates of total richness (N) and estimates of richness at each of the 
#J sampling locations (Nsite)
N = fit$sims.list$N
mean(N); summary(N); plot(table(N))

Nsite = fit$sims.list$Nsite
site.richness.matrix = cbind(apply(Nsite,2,mean), apply(Nsite,2,mean))

#Plot mean site richness against one of the covariates to examine how point 
#richness varies as a result of understory foliage
plot(ufc, apply(Nsite,2,mean), pch=16, lwd=2, xlab="Understory foliage (UFC)",
        ylab="Point richness", type="p")

#Examine how mean species-specific occupancy changes by ufc in one of the 
#habitat types (CAT0)
x = seq(-1.5,2.5, by=0.1)
y = (x*sdufc) + mufc
a1 = fit$sims.list$a1 
a2 = fit$sims.list$a2

plot(y,y, type="l", ylim=c(0,1), xlim=c(0,40), main="Species-specific relationships with ufc",
col="white")

for (i in 1:n) {
   CATO <- mean(cato.occ[,i]) + mean(a1[,i])*x + mean(a2[,i])*x*x
   lines(y,plogis(CATO), type="l", ylim=c(0,1), xlim=c(0,40),   
		main=i)
   }
