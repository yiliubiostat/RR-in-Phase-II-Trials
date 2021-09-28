# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~ BIOSTAT 907 Term Project R Code ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~ Function for Simon's Optimal and Minimax Design ~~~~~
# ~~~~~~~~ with Futility Stopping Only                     ~~~~~

# Author: Yi Liu
# Create date: 09/23/2021

# This script provides a function of calculating 
# --- Simon's optimal two stage design (1989): to minimize the expected sample size EN
# --- and minimax two stage design: to minimize the "maximal" sample size n1 + n2

# The designs in this script are with only the futility stopping (FS)


# Working path using Duke Computer Cluster (DCC)
# Don't run this on your server; please modify the path if necessary
# setwd("/hpc/home/yl724")


design_FS <- function(desired.alpha, desired.power, 
                   p0, p1, max.N = 55, type = "optimal") {
  
  # Testing: H0: p = p0 vs. H1: p = p1
  # alpha: the type I error rate
  # --- typically, alpha = 0.1, 0.05, etc.
  # power: = 1 - beta, beta is the type II error rate
  # --- typically, power = 0.8, 0.85, 0.9, etc.
  # p0: the response rate of the standard (historical) therapy
  # p1: the response rate of the experimental therapy
  # --- typically, p1 - p0 = 0.2
  # max.N: the maximal tolerance of feasible sample size, is 55 by default from experience
  # type: can be "optimal" or "minimax", by default is "optimal"
  
  
  if(p0 < 0 | p0 > 1 | p1 < 0 | p1 > 1) stop("p0 and p1 must be in the range [0, 1]!")
  if(p1 <= p0) stop("p1 must be greater than p0; we only consider one-sided test!")
  
  if(max.N <= 15) stop("The smallest sample size must be greater than 16!")
  
  if(max.N > 100) stop("Tha maximal sample size is too large for a phase II trial, we can't do this!")
  
  if(!(type %in% c("optimal", "minimax"))) stop("Design can only either be optimal or minimax!")
  
  # Initial values
  
  EN0 = max.N
  
  S1.size = 0
  S1.bound = 0
  
  total.size = 0
  S2.bound = 0
  
  alpha = NA
  power = NA
  
  # -------------------------------------------------------------------------------
  # Find minimax design

  if(type == "minimax") {
    
    for(n in 16:max.N) {
      for(n1 in 1:(n-1)) {
        for(a1 in 0:n1) {
          for(a in (a1+1):n) {
            
            if(a1+1 <= n1) {
              
              type.I.err <- sum(
                dbinom(x = (a1+1):n1, size = n1, prob = p0)*
                  (1-pbinom(q = (a-a1-1):(a-n1), size = n-n1, prob = p0)) )
              
              stat.pwr <- sum(
                dbinom(x = (a1+1):n1, size = n1, prob = p1)* 
                  (1-pbinom(q = (a-a1-1):(a-n1), size = n-n1, prob = p1)) )
              
              PET = pbinom(q = a1, size = n1, prob = p0)
              EN = PET*n1 + (1-PET)*n
            
            } else { next }
            
            if(type.I.err <= desired.alpha & stat.pwr >= desired.power & EN < EN0) {
              
              alpha = type.I.err
              power = stat.pwr
              
              S1.bound = a1
              S1.size = n1
              
              S2.bound = a
              total.size = n
              
              EN0 = EN
              
            } else { next }
          }
        }
      }
      if(is.na(alpha) == FALSE) { break }
      n = n+1
    }
    
    if(is.na(power)) { warning("No design found! Cannot obtain enough power by assigned parameters and conditions!") }
    
    minimax.design <- data.frame(p0 = p0, p1 = p1,
                                 a1 = S1.bound, N1 = S1.size,
                                 a = S2.bound, N = total.size,
                                 EN = round(EN0, 2), 
                                 alpha = alpha, power = power)
    return(minimax.design)
  }
  
  
  # ------------------------------------------------------------------------------
  # Find optimal design
  
  if(type == "optimal") {
    
    for(n in 16:max.N) {
      for(n1 in 1:(n-1)) {
        for(a1 in 0:n1) {
          for(a in (a1+1):n) {
            
            if(a1+1 <= n1) {
              
              type.I.err <- sum(
                dbinom(x = (a1+1):n1, size = n1, prob = p0)*
                  (1-pbinom(q = (a-a1-1):(a-n1), size = n-n1, prob = p0)) )
              
              stat.pwr <- sum(
                dbinom(x = (a1+1):n1, size = n1, prob = p1)* 
                  (1-pbinom(q = (a-a1-1):(a-n1), size = n-n1, prob = p1)) )
              
              PET = pbinom(q = a1, size = n1, prob = p0)
              EN = PET*n1 + (1-PET)*n
              
            } else { next }
            
            if(type.I.err <= desired.alpha & stat.pwr >= desired.power & EN < EN0) {
              
              alpha = type.I.err
              power = stat.pwr
              
              S1.bound = a1
              S1.size = n1
              
              S2.bound = a
              total.size = n
              
              EN0 = EN
              
            } else { next }
          }
        }
      }
      n = n+1
    }
    
    if(is.na(power)) { warning("No design found! The difference of p1 and p0 is so small to obtain enough power.") }
    
    optimal.design <- data.frame(p0 = p0, p1 = p1,
                                 a1 = S1.bound, N1 = S1.size,
                                 a = S2.bound, N = total.size,
                                 EN = round(EN0, 2), 
                                 alpha = alpha, power = power)
    return(optimal.design)
  }
  
}



# ~~~~~~~~~~~~~~~~~~~~~~ END of the program ~~~~~~~~~~~~~~~~~~~~~~~~~



