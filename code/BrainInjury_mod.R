
## The data in this folder was provided by 
# Genest 
# It comprises processed data, including data used in:
# 
#   C Genest, J G Nešlehová, B Rémillard, O A Murphy, 
# Testing for independence in arbitrary distributions, Biometrika, Volume 106, Issue 1, March 2019, 
# Pages 47–68, https://doi.org/10.1093/biomet/asy059

BrainInjury <- read.csv("/Volumes/GoogleDrive/Il mio Drive/autodep/R/Code/Data/BrainInjury.csv",sep=",",na=c("n/a"))

attach(BrainInjury)
# Checking for missing values 
which(is.na(Education))
# [1] 5
which(is.na(BNT))
# [1] 139
which(is.na(FASmean) & is.na(P))
# [1] 10 25
which(is.na(vetements) & is.na(animals))
which(is.na(BDAE))
# [1] 105
which(is.na(SCATBIcor))
which(is.na(DTLA))
which(is.na(ChapmanCorrect))
# [1]  5 10 14 17 28 45 46

# Removing missing values of the variables we are interested in (Age, Education, GCS,LOS, GOSE, BNT, DTLA, MEC,Litteral fluency, categorial fluency,BDAE,SCATBIcor)

Ind <- rep("TRUE",145)
Ind[c(5,10,14,17,25,28,45,46,105,139)] <- FALSE

#detach(BrainInjury)
B <- asNumericMatrix(BrainInjury)
BS <- B[c(1:4,6:9,11:13,15:16,18:24,26:27,29:44,47:104,106:138,140:145),c(1:35)]
colnames(BS)

# [1] "MRN"                 "Sexe"                "Age"                 "Education"           "Langue"              "LOS"                
# [7] "Assurance"           "Transfert"           "GCS"                 "Severite"            "Perte.de.conscience" "GOSE"               
# [13] "Marshall"            "Frontal.G"           "TemporalG"           "ParietalG"           "FrontalD"            "TemporalD"          
# [19] "ParietalD"           "Occipital"           "DTLA"                "BNT"                 "ChapmanCorrect"      "ChapmanError"       
# [25] "ChapmanTotal"        "MECconv"             "P"                   "vetements"           "BDAE"                "SCATBIbrut"         
# [31] "SCATBIcor"           "TRFB"                "TRFC"                "animals"             "FASmean"   


#Pooling together "FASmean" and "P"

LF <- as.numeric(135)
for(i in 1:135){
  if(BS[i,5]==2){LF[i] <- BS[i,35]}
  if(BS[i,5]==1){LF[i] <- BS[i,27]}
}

#Pooling together "animals" and "vetements"


CF <- as.numeric(135)
for(i in 1:135){
  if(BS[i,5]==2){CF[i] <- BS[i,34]}
  if(BS[i,5]==1){CF[i] <- BS[i,28]}
}

# Binning GOSE

GOSEBin <- as.numeric(BS[,12]<5)   

# Only 19 patients with severe GOSE! GOSE then appears independent of the language tests


# Bn test  ----------------------------------------------------------------

setwd("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2")

#source("indeptest.R")
source("indeptest_with_cat.R")
source("libraries.R")
source("utilities.R")

mydata <- cbind(BS[,c(3,4,12,21,22,23,26,29,9)],LF,CF)
dim(mydata)

str(mydata)
varnames <- c(colnames(BS[,c(3,4,12,21,22,23,26,29,9)]),"LF","CF")
varnames[2] <- "Edu"
varnames[7] <- "MEC"
varnames[6] <- "CC"

colnames(mydata) <- varnames

head(mydata)

# Age <- mydata[,1] ## Age
# Edu <- mydata[,2] ## Edu
# GCS <- mydata[,9] ## GCS
# GOSE <- mydata[,3] ## GOSE

## language tests 
# ltest <- mydata[,-c(1,2,3,9)]
# dim(ltest)

p.order <- function(n) max(1, floor(n^(1/3))-1)

Bn.pval <- matrix(NA, ncol=11, nrow=11)
#Bn.pval.p1 <- matrix(NA, ncol=11, nrow=11)
for(i in 1:11){
  for(j in 1:11){
    x <- as.integer(mydata[,i])
    y <- as.integer(mydata[,j])
    test <- try(indeptest(x,y, ties.method = "first")$pvalue, silent = T)
    #test <- try(indeptest(x,y)$pvalue, silent = T)
    Bn.pval[i,j]    <- ifelse(is.numeric(test), test, NA)
  #   Bn.pval.p1[i,j] <- indeptest(mydata[,i],mydata[,j],
  #                                order=c(p.order,1), ties.method = "first")$pvalue
    }
}


colnames(Bn.pval) <- rownames(Bn.pval) <- varnames
round(Bn.pval,3)

index<- upper.tri(Bn.pval, diag = FALSE)
Bn.pval[lower.tri(Bn.pval)] <- NA
Bn.pval

Bn.pval_BH <- round(p.adjust(na.omit(as.vector(Bn.pval)), 
                             method = "BH"),3)
length(Bn.pval_BH)

Bn.pval_BH.matrix <- matrix(NA, nrow=11, ncol=11,
                            byrow=T)

Bn.pval_BH.matrix[1, ] <- c(1,Bn.pval_BH[1:10]) # 10
Bn.pval_BH.matrix[2, ] <- c(NA,1,Bn.pval_BH[11:19]) # 9
Bn.pval_BH.matrix[3, ] <- c(NA, NA,1,Bn.pval_BH[20:27])#8
Bn.pval_BH.matrix[4, ] <- c(NA, NA, NA,1,Bn.pval_BH[28:34])# 7
Bn.pval_BH.matrix[5, ] <- c(NA, NA, NA, NA,1,Bn.pval_BH[35:40])#6
Bn.pval_BH.matrix[6, ] <- c(NA, NA, NA, NA, NA,1,Bn.pval_BH[41:45])#5
Bn.pval_BH.matrix[7, ] <- c(NA,NA, NA, NA, NA, NA,1,Bn.pval_BH[46:49])#4
Bn.pval_BH.matrix[8, ] <- c(NA, NA, NA, NA, NA, NA, NA,1,Bn.pval_BH[50:52])#3
Bn.pval_BH.matrix[9, ] <- c(NA, NA, NA, NA, NA, NA, NA, NA,1,Bn.pval_BH[53:54])#2
Bn.pval_BH.matrix[10, ] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA,1,Bn.pval_BH[55])#1
Bn.pval_BH.matrix[11, ] <- c(rep(NA,10),1)

Bn.pval_BH.matrix[lower.tri(Bn.pval_BH.matrix)]<-
  t(Bn.pval_BH.matrix)[lower.tri(Bn.pval_BH.matrix)]

colnames(Bn.pval_BH.matrix) <- rownames(Bn.pval_BH.matrix) <- varnames
Bn.pval_BH.matrix



library(xtable) 
xtable(Bn.pval_BH.matrix)


# "Age"  "Edu"  "GOSE" "DTLA" "BNT" 
# "CC"   "MEC"  "BDAE" "GCS"  "LF"   "CF" 

sum(10:1)
# [1] 55  # choose(11,2)

# explore pairwise dependence
# with the Benjamini–Hochberg adjustment for false discovery rate;
# pairwise independence that were significant at the 5% level.

FDR=.05
significativi <- function(x) sum(x<=FDR)

apply(Bn.pval_BH.matrix,1,significativi)

# "Age"  "Edu"  "GOSE" "DTLA" "BNT" 
# "CC"   "MEC"  "BDAE" "GCS"  "LF"   "CF" 

# Age
which(Bn.pval_BH.matrix[1, 2:11]<=FDR) # BNT  MEC BDAE  GCS   CF  
# Edu
which(Bn.pval_BH.matrix[2, 3:11]<=FDR) # BNT  CC MEC  LF  CF 
# GOSE
which(Bn.pval_BH.matrix[3, 4:11]<=FDR) # DTLA  BNT BDAE  GCS   LF
# DTLA
which(Bn.pval_BH.matrix[4, 5:11]<=FDR) # BNT   CC BDAE  GCS 
# BNT
which(Bn.pval_BH.matrix[5, 6:11]<=FDR) # CC MEC  CF  
# CC
which(Bn.pval_BH.matrix[6, 7:11]<=FDR) # MEC BDAE  GCS 
# MEC
which(Bn.pval_BH.matrix[7, 8:11]<=FDR) # CF 
# BDAE
which(Bn.pval_BH.matrix[8, 9:11]<=FDR) # GCS  LF  CF 
# "GCS"  
which(Bn.pval_BH.matrix[9, 10:11]<=FDR)
# "LF"  
which(Bn.pval_BH.matrix[10, 11]<=FDR)


# Age, 
# education, Edu, 
# severity of injury, GCS,
# level of disability at discharge, GOSE, 
# seven language tests: 
#   the Boston Diagnostic Aphasia Examination, BDAE, 
#   the Boston Naming Test, BNT, 
#   the Chapman–Cook Speed of Reading Test, CC, 
#   the Categorial Fluency test, CF, 
#   the Detroit Test of Learning Aptitude, DTLA, 
#   the Protocole Montréal d’évaluation de la communication, MEC, 
#   and the Literal Fluency test, LF.





