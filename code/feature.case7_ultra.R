
rm(list=ls())

setwd("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2")

source("indeptest_with_cat.R")
source("libraries.R")
source("utilities.R")
source("gendataRegMod.R")


type.order<- function(n){
  # max(2, floor((n)^(1/2)) - 3)
  max(1, floor(n^(1/3))-1)
}


nSim     <- 1000
p        <- 2500
nsample  <- 100
set.seed(123)


### Selection Mats
pMat = numeric(nSim)
dMat = numeric(nSim)
BnMat = numeric(nSim)
BnMat_pq = numeric(nSim)


active <- 1:4

for(i in 1:nSim){
  ### Generate Y and X
  sim.data <- Gendata2f(n=nsample, p=p, rho=0.5)
  X0 <- sim.data$X
  Y0 <- sim.data$Y
  
  ########### Pearson Correlation ##########
  pCorrs = abs(cor(X0,Y0,method = "pearson"))
  pCorrOrder = order(pCorrs, decreasing = T)
  pMat[i] <- max(which(pCorrOrder %in% active))
  
  ########### Distance Correlation #########
  dCorrs = abs(apply(X0, 2, dcor, y = Y0))
  dCorrOrder = order(dCorrs, decreasing = T)
  dMat[i] <- max(which(dCorrOrder %in% active))
  
  ########### Bn-Screen p1  #############
  indepy <- function(Xk) {
    indeptest(Xk, Y0, order =c(type.order(nsample),1), basis = "spline")$stat
  }
  BnCorrs <- abs(apply(X0, 2, indepy))
  BnCorrOrder = order(BnCorrs, decreasing = T)
  BnMat[i] <- max(which(BnCorrOrder %in% active))
  
  ########### Bn-Screen  pq #############
  indepy_pq <- function(Xk) {
    indeptest(Xk, Y0, basis = "spline")$stat
  }
  BnCorrs_pq <- abs(apply(X0, 2, indepy_pq))
  BnCorrOrder_pq = order(BnCorrs_pq, decreasing = T)
  BnMat_pq[i] <- max(which(BnCorrOrder_pq %in% active))
  
  print(i)
}



## Average screening proportions of true signals

## The quantiles of minimum model size for linear and Poisson models in Example 1
# over 1000 replications. The true model size is 5.

out.1a <- cbind(pMat, dMat, BnMat, BnMat_pq)

nome <- "case_7_ultra"
save(out.1a,
     file = paste0("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/results/varsel/", nome, ".Rdata"))

# res.1a <- cbind(apply(out.1a, 2, quantile, probs  = 0.05),
# apply(out.1a, 2, quantile, probs  = 0.25),
# apply(out.1a, 2, quantile, probs  = 0.5),
# apply(out.1a, 2, quantile, probs  = 0.75),
# apply(out.1a, 2, quantile, probs  = 0.95))
# 
# colnames(res.1a) <- c("5%","25%", "50%", "75%", "95%")
# res.1a

#boxplot(pMat, dMat, BnMat)
# 
# 
# 

