
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
      
