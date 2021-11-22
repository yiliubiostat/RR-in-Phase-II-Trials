# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~ BIOSTAT 907 Term Project R Code ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~ Trends of RR Boundary with Futility Stopping ~~~~~~~~~

# Author: Yi Liu
# Create date: 09/25/2021

# Working path using Duke Computer Cluster (DCC)
# Don't run this on your server; please modify the path if necessary
# setwd("/hpc/home/yl724")

source("design_FS_func.R")

data_FS <- function(alpha, power, diff=0.2, p0.min=0.05, p0.max=0.7, p0.incre=0.005) {
  
  # alpha: desired type I error level for each design
  # power: desired power for each design
  # diff: = p1-p0, typically assigned to be 0.2 ; suggested >= 0.15
  
  P0 <- seq(p0.min, p0.max, by=p0.incre)
  P1 <- P0 + diff
  P <- (P0 + P1)/2
  l <- length(P0)
  
  op.design <- data.frame(p0 = NA, p1 = NA,
                          a1 = NA, N1 = NA, a = NA, N = NA, EN = NA,
                          alpha = NA, power = NA)
  mm.design <- data.frame(p0 = NA, p1 = NA,
                          a1 = NA, N1 = NA, a = NA, N = NA, EN = NA,
                          alpha = NA, power = NA)
  
  for(p0 in P0) {
    p1 = p0 + diff
    
    # call the function defined before, and set the maximal sample size to be 70
    op.design <- rbind(op.design, design_FS(alpha, power, p0 = p0, p1 = p1, type = "optimal", max.N = 55))
    mm.design <- rbind(mm.design, design_FS(alpha, power, p0 = p0, p1 = p1, type = "minimax", max.N = 55))
    
  }
  return(list(P0 = P0, P1 = P1, P = P, 
              op.design = op.design[-1,], 
              mm.design = mm.design[-1,]))
}


# We generate data based for choosing 
# --- (alpha, power) = (0.1, 0.8), (0.1, 0.85), (0.1, 0.9), (0.05, 0.8), (0.05, 0.85), (0.05, 0.9)
# --- diff = 0.2, 0.25

ls1.1 <- data_FS(alpha = 0.1, power = 0.8, diff = 0.2)
ls1.2 <- data_FS(alpha = 0.1, power = 0.8, diff = 0.25)

ls2.1 <- data_FS(alpha = 0.1, power = 0.85, diff = 0.2)
ls2.2 <- data_FS(alpha = 0.1, power = 0.85, diff = 0.25)

ls3.1 <- data_FS(alpha = 0.1, power = 0.9, diff = 0.2)
ls3.2 <- data_FS(alpha = 0.1, power = 0.9, diff = 0.25)

ls4.1 <- data_FS(alpha = 0.05, power = 0.8, diff = 0.2)
ls4.2 <- data_FS(alpha = 0.05, power = 0.8, diff = 0.25)

ls5.1 <- data_FS(alpha = 0.05, power = 0.85, diff = 0.2)
ls5.2 <- data_FS(alpha = 0.05, power = 0.85, diff = 0.25)

ls6.1 <- data_FS(alpha = 0.05, power = 0.9, diff = 0.2)
ls6.2 <- data_FS(alpha = 0.05, power = 0.9, diff = 0.25)


save(file = "design_FS.Rdata", 
     ls1.1, ls1.2, ls2.1, ls2.2,
     ls3.1, ls3.2, ls4.1, ls4.2,
     ls5.1, ls5.2, ls6.1, ls6.2)


# ~~~~~~~~~~~~~~~~~~~~~~ END of the program ~~~~~~~~~~~~~~~~~~~~~~~~~


