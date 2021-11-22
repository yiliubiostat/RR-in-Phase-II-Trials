# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~ BIOSTAT 907 Term Project R Code ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~ Function for Simon's Optimal and Minimax Design ~~~~~
# ~~~~~~~~ with both Futility and Superiority Stopping     ~~~~~

# Author: Yi Liu
# Create date: 10/11/2021

# This script provides a function of calculating 
# --- Simon's optimal two stage design (1989): to minimize the expected sample size EN
# --- and minimax two stage design: to minimize the "maximal" sample size n1 + n2

# The designs in this script are with both futility and superiority stopping


# Working path using Duke Computer Cluster (DCC)
# Don't run this on your server; please modify the path if necessary
# setwd("/hpc/home/yl724")


design_twoS <- function(desired.alpha, desired.power, 
                      p0, p1, max.N = 55, type = "optimal") {
  
  # Testing: H0: p = p0 vs. H1: p = p1 (Simple vs. Simple)
  # alpha: the type I error rate
  # --- typically, alpha = 0.1, 0.05, etc.
  # power: = 1 - beta, beta is the type II error rate
  # --- typically, power = 0.8, 0.85, 0.9, etc.
  # p0: the response rate of the standard (historical) therapy
  # p1: the response rate of the experimental therapy
  # --- typically, p1 - p0 = 0.2
  # max.N: the maximal tolerance of feasible sample size, is 55 by default from experience
  # type: type of design; can be "optimal" or "minimax", by default is "optimal"
  
  
  if(p0 < 0 | p0 > 1 | p1 < 0 | p1 > 1) stop("p0 and p1 must be in the range [0, 1]!")
  if(p1 <= p0) stop("p1 must be greater than p0; we only consider one-sided test!")
  
  if(max.N <= 15) stop("The smallest sample size must be greater than 16!")
  
  if(max.N > 100) stop("Tha maximal sample size is too large for a phase II trial, we can't do this!")
  
  if(!(type %in% c("optimal", "minimax"))) stop("Design can only either be optimal or minimax!")
  
  # Initial values
  
  EN.init = max.N
  
  S1.size = 0 # n1
  S1.rej = 0 # a1
  S1.cont = 0 # b1
  
  total.size = 0 # n (n1+n2 =n )
  S2.bound = 0 # a
  
  alpha = NA
  power = NA
  
  # -------------------------------------------------------------------------------
  # Find minimax design
  
  if(type == "minimax") {
    
    for(n in 16:max.N) {
      for(n1 in 1:(n-1)) {
        for(a1 in 0:n1) {
          for(b1 in (a1+1):n1) {
            for(a in b1:n) {
              
              if(a1+1 <= b1-1) {
                
                tp.I.err <- 1 - (pbinom(q = a1, size = n1, prob = p0) + 
                                   sum(dbinom(x = (a1+1):(b1-1), size = n1, prob = p0)*
                                         pbinom(q = (a-a1-1):(a-b1+1), size = n-n1, prob = p0)) )
                
                stat.pwr <- 1 - (pbinom(q = a1, size = n1, prob = p1) + 
                                   sum(dbinom(x = (a1+1):(b1-1), size = n1, prob = p1)*
                                         pbinom(q = (a-a1-1):(a-b1+1), size = n-n1, prob = p1)) )
                
                PET0 = pbinom(q = a1, size = n1, prob = p0) + 1 - pbinom(q = b1-1, size = n1, prob = p0)
                PET1 = pbinom(q = a1, size = n1, prob = p1) + 1 - pbinom(q = b1-1, size = n1, prob = p1)
                
                EN0 = PET0*n1 + (1-PET0)*n
                EN1 = PET1*n1 + (1-PET1)*n
                EN = mean(c(EN0, EN1))
                
              } else { next }
              
              if(tp.I.err <= desired.alpha & stat.pwr >= desired.power & EN < EN.init) {
                
                alpha = tp.I.err
                power = stat.pwr
                
                S1.rej = a1
                S1.cont = b1
                S1.size = n1
                
                S2.bound = a
                total.size = n
                
                EN.init = EN
                
              } else { next }
            }
          }
        }
      }
      if(is.na(alpha) == FALSE) { break }
    }
    
    if(is.na(power)) { warning("No design found! Cannot obtain enough power by assigned parameters and conditions!") }
    
    minimax.design <- data.frame(p0 = p0, p1 = p1,
                                 a1 = S1.rej, b1 = S1.cont, N1 = S1.size,
                                 a = S2.bound, b = S2.bound+1, N = total.size,
                                 EN = round(EN.init, 2), 
                                 alpha = alpha, power = power)
    return(minimax.design)
  }
  
  
  # ------------------------------------------------------------------------------
  # Find optimal design
  
  if(type == "optimal") {
    
    for(n in 16:max.N) {
      for(n1 in 1:(n-1)) {
        for(a1 in 0:n1) {
          for(b1 in (a1+1):n1) {
            for(a in b1:n) {
              
              if(a1+1 <= b1-1) {
                
                tp.I.err <- 1 - (pbinom(q = a1, size = n1, prob = p0) + 
                                   sum(dbinom(x = (a1+1):(b1-1), size = n1, prob = p0)*
                                         pbinom(q = (a-a1-1):(a-b1+1), size = n-n1, prob = p0)) )
                
                stat.pwr <- 1 - (pbinom(q = a1, size = n1, prob = p1) + 
                                   sum(dbinom(x = (a1+1):(b1-1), size = n1, prob = p1)*
                                         pbinom(q = (a-a1-1):(a-b1+1), size = n-n1, prob = p1)) )
                
                PET0 = pbinom(q = a1, size = n1, prob = p0) + 1 - pbinom(q = b1-1, size = n1, prob = p0)
                PET1 = pbinom(q = a1, size = n1, prob = p1) + 1 - pbinom(q = b1-1, size = n1, prob = p1)
                
                EN0 = PET0*n1 + (1-PET0)*n
                EN1 = PET1*n1 + (1-PET1)*n
                EN = mean(c(EN0, EN1))
                
              } else { next }
              
              if(tp.I.err <= desired.alpha & stat.pwr >= desired.power & EN < EN.init) {
                
                alpha = tp.I.err
                power = stat.pwr
                
                S1.rej = a1
                S1.cont = b1
                S1.size = n1
                
                S2.bound = a
                total.size = n
                
                EN.init = EN
              
              }  else { next }
            }
          }
        }
      }
    }
    
    if(is.na(power)) { warning("No design found! Cannot obtain enough power by assigned parameters and conditions!") }
    
    optimal.design <- data.frame(p0 = p0, p1 = p1,
                                 a1 = S1.rej, b1 = S1.cont, N1 = S1.size,
                                 a = S2.bound, b = S2.bound+1, N = total.size,
                                 EN = round(EN.init, 2), 
                                 alpha = alpha, power = power)
    return(optimal.design)
  }
  
}



# ~~~~~~~~~~~~~~~~~~~~~~ END of the program ~~~~~~~~~~~~~~~~~~~~~~~~~



