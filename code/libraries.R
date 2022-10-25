
##############################################################################
###
### required libraries for comparison  
###
##############################################################################

library(copula)           # for indepTest (Genest and Rémillard, 2004)
library(energy)           # for dCov (Szekely et al., 2007)
library(geigen)           # for generalized eigenvalues used in indeptest
library(ggplot2)
library(gridExtra)
library(HHG)              # for HHG test
library(Hmisc)            # for Hoeffding's test
library(microbenchmark)   # to test for speed
library(minerva)          # for the MIC
library(splines)          # for the spline basis
