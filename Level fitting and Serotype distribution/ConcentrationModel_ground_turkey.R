#####################################################################
## Modeling MPN concentrations in JAGS (called from R using rjags)

## Authors: Elisabetta Lambertini (RTI International). Last edited Aug 2019.
## Adapted from OpenBUGS code from Mike Williams, FSIS (Williams and Ebel, 2012).
## This scripts models all samples that were screened (presence/absence).
#####################################################################

# Empty global environment and load libraries
rm(list = ls())
library(lattice)
library(rjags)
library(pastecs)
library(coda)

getwd()
# set working directory:
setwd("C:/Users/yiyili2/Box/Poultry/USPoultry work part 2/data and analysis/Fitting")

# Assumptions:
# 1) each sample analyzed is a realization of a Poisson distribution (not overinflated)
# 2) All samples in a set are considered here: those negative (non-detected) at screening, those positive at screening and not further anayzed, and those positive for screening and then enumerated via MPN.
# 3) There are no "true negative" samples. Samples negative at screening are considered to have an average concentration that makes them unlikely to be detected.

# Specify data:

# Parameters for sampling results (screening and MPN assay):
N.total = 4284       # overall number of samples tested
N.pos.test = 766     # number of screening-positive (detected) samples. Not all were enumerated.
N.pos.test.mpn = 179 # number of samples positive (detected) at +/- screening test that were also enumerated via MPN
num.tubes = 3        # number of tubes at each dilution

#tube2 <-read.csv("Tube_scores.csv", head=F) # if loading a data file

# Data:  # These are MPN tube scores from the combined "NRTE" and "HC_TU" datasets
tube2 <- matrix(c(
  0,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,2,0,0,0,
  1,0,0,0,0,
  3,2,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  3,3,0,1,0,
  0,0,0,0,0,
  1,0,0,0,0,
  1,1,0,0,0,
  1,0,0,0,0,
  1,2,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  3,1,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  1,1,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  1,1,0,0,0,
  3,2,1,0,0,
  2,0,0,0,0,
  1,0,0,0,0,
  3,1,0,0,0,
  0,0,0,0,0,
  1,0,1,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  2,0,0,0,0,
  3,0,0,0,0,
  2,1,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  3,3,0,0,0,
  3,3,0,1,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,3,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  2,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  0,0,0,0,0,
  3,1,0,0,0,
  2,0,0,0,0,
  0,0,0,0,0,
  3,3,0,1,0,
  1,0,0,0,0,
  3,0,0,0,0,
  3,2,0,0,0,
  0,0,0,0,0,
  3,3,1,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,0,0,0,0,
  3,1,0,0,0,
  1,3,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,2,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,1,0,0,0,
  1,0,0,0,0,
  3,3,1,0,0,
  3,3,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,0,0,0,0,
  0,0,0,0,0,
  2,0,1,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  3,0,0,0,0,
  2,0,0,0,0,
  3,0,0,0,0,
  2,1,0,0,1,
  2,2,0,0,0,
  3,0,0,0,0,
  2,0,0,0,0,
  0,0,0,0,0,
  1,1,0,0,0,
  1,0,0,0,0,
  3,0,0,0,0,
  3,3,1,0,0,
  3,3,3,1,0,
  0,0,0,0,0,
  3,3,1,0,0,
  2,0,0,0,0,
  3,1,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,2,0,0,0,
  3,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,1,0,0,0,
  2,1,0,0,0,
  0,0,0,0,0,
  3,2,0,0,0,
  1,0,0,0,0,
  1,0,0,0,0,
  1,0,0,0,0,
  1,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  3,2,1,1,0,
  1,0,0,0,0,
  1,0,0,0,0,
  3,3,2,0,0,
  0,0,0,0,0,
  3,3,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  2,0,0,0,0,
  0,0,0,0,0,
  3,3,3,3,0,
  0,0,0,0,0,
  0,0,0,0,0),nrow=179,ncol=5,byrow=T)

v.screen = 325 # screen test volume (ml=g)
poscr     <- rep(1,N.pos.test.mpn) # add column for screening results (here assumed all to be 1, as only screen-positives were enumerated via MPN)
tube2plus <- cbind(poscr,tube2)

screen.results <- c(rep(0,N.total-N.pos.test),rep(1,N.pos.test-N.pos.test.mpn))

# JAGS model:

{ # Extra bracket needed only for R markdown files
  sink("MPN_JAGSmodel_v4.R")
  cat("
      model{
      
      # PRIORS
      # Broad priors of a lognormal describing CFU per gram of product.  
      # (Note: if using Win/OpenBugs, they use the natural log (ln) scale so also include a transformation to the log base 10 scale).
      
      mu ~ dnorm(-9.4,1.1)
      tau ~ dgamma(0.1,0.1)
      sigma <- sqrt(1 / tau)
      
      # Loop for the screen-positive tests that were also evaluated by MPN:
      
      for (i in 1:N.pos.test.mpn) {
      
      lambda[i] ~ dlnorm(mu,tau)
      p_screen[i] <- 1- exp(-v.screen*lambda[i]) # 325 screening test
      p_10[i]     <- 1- exp(-10.0*lambda[i])     # Probabilities for MPN tubes
      p_1[i]      <- 1- exp(-1.0*lambda[i])
      p_0.1[i]    <- 1- exp(-0.1*lambda[i])
      p_0.01[i]   <- 1- exp(-0.01*lambda[i])
      p_0.001[i]  <- 1- exp(-0.001*lambda[i])
      
      # Likelihood for the qualitative screening test
      tube[i,1] ~ dbin(p_screen[i],1) # Test positive on 325 g screening test
      
      # Likelihoods for the MPN data at each dilution
      tube[i,2] ~ dbin(p_10[i],N.tubes)  # with N.tunes=3, the draw is a number btw 0 and 3, included.
      tube[i,3] ~ dbin(p_1[i],N.tubes)
      tube[i,4] ~ dbin(p_0.1[i],N.tubes)
      tube[i,5] ~ dbin(p_0.01[i],N.tubes)
      tube[i,6] ~ dbin(p_0.001[i],N.tubes)
      
      }  # end of loop on screen-positive MPN-enumerated samples (i)
      
      
      # Loop for the screened samples that had no MPN performed 
      # (both detected and non-detected at +/- screening test):
      
      for (j in (N.pos.test.mpn+1):N.total){   
      lambda[j] ~ dlnorm(mu,tau)
      p_screen[j] <- 1-exp(-v.screen*lambda[j])
      screened[j-N.pos.test.mpn] ~ dbin(p_screen[j],1)
      } 
      
      
      
      }  # end of model
      ", fill = TRUE)
  sink()
} # Extra bracket needed only for R markdown files


# Specify data and inits:
data = list(
  v.screen = v.screen,
  N.tubes = num.tubes,
  N.total = N.total,
  N.pos.test.mpn = N.pos.test.mpn,
  tube = tube2plus,
  screened = screen.results
) 

inits = list(
  list(lambda = rlnorm(N.total,-3.9,0.3),tau = 1.1,.RNG.seed=1,.RNG.name= 'base::Wichmann-Hill'), 
  list(lambda = rlnorm(N.total,-3.0,0.1),tau = 1.1,.RNG.seed=2,.RNG.name= 'base::Wichmann-Hill'), 
  list(lambda = rlnorm(N.total,-4,0.1),  tau = 1.1,.RNG.seed=3,.RNG.name= 'base::Wichmann-Hill') 
)

# Specification for the MCMC sampler in JAGS:
n.adapt   = 1000   # n iteration for adaptation
n.update  = 10000  # n iterations for burn-in
n.iter    = 50000  # n iteration to keep in the final chain
n.thinval = 50     # Thinning constant (50).  (Use n.thinval = 1 when just checking out the model run).
burn.in   = 1000   # 1000

# Call to JAGS:
set.seed(1)
jm_a = jags.model("MPN_JAGSmodel_v4.R", n.adapt = n.adapt, n.chains = 3, inits = inits, data = data)
#update(jm_a, n.iter = n.update)  # repeat this line if more iteration are needed for convergence
zm1 = coda.samples(jm_a, variable.names = c("mu","tau"), n.iter = n.iter, thin=n.thinval) 

(SummaryResults=summary(window(zm1, start = burn.in), quantiles=c(0.025, 0.5, 0.975)))
mu=    c(window(zm1[[1]][,1],start=burn.in),
         window(zm1[[2]][,1],start=burn.in),
         window(zm1[[3]][,1],start=burn.in))

tau=   c(window(zm1[[1]][,2],start=burn.in),
         window(zm1[[2]][,2],start=burn.in),
         window(zm1[[3]][,2],start=burn.in))
sigma= sqrt(1/tau)

# Summarizing coda object:
summary(zm1)  # print(zm1)  to visualize the entire output
summary(tau)
#plot(zm1)
traceplot(zm1)

# Test for convergence:
heidel.diag(zm1)

### Plots and outputs:
summary_stats<-as.data.frame(cbind(stat.desc(mu),stat.desc(sigma)))  
names(summary_stats)<-c("mu","sigma")
write.csv(summary_stats,"JAGS_Fit_Summary.csv")

# Plot density functions:
mu<-as.matrix(sort(mu))
K.ecdf<-ecdf(mu)
par(mfrow=c(1,1),mar=c(3,3,2,2), las=1)
plot(density(mu), main="mu")

sigma<-as.matrix(sort(sigma))
K.ecdf<-ecdf(sigma)
plot(density(sigma), main="sigma")

par(mfrow=c(1,2),mar=c(1,3,2,2), las=1)
boxplot(mu,main="mu",col="lightgreen")
boxplot(sigma,main="sigma",col="yellow")

# Print summary statistics of fitted parameters of the lognormal distribution:
print(summary_stats)

###################


