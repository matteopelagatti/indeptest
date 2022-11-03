# Cardiomyopathy microarray data ------------------------------------------------------------------

rm(list=ls())

setwd("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2")

#source("indeptest.R")

library(gamsel)

source("indeptest_with_cat.R")
source("libraries.R")
source("utilities.R")

# These data were used by Segal, Dahlquist, and Conklin (2003) to evaluate regression-based approaches to microarray
# analysis. The aim was to determine which genes were influential for overexpression of
# a G protein-coupled receptor, designated Rol, in mice. The research related to understand ing types of human heart disease. The Rol expression level, F/, was measured for n = 30
# specimens, and genetic expression levels, X;, were obtained for p = 6,319 

# Read in data
# X <- t(read.table("/Users/giannamonti/Downloads/hall-miller-supplements/Data/Ro131.csv", header=TRUE, sep=","))
# dim(X)

X <- t(read.table("Data/Ro131.csv",
                  header=TRUE, sep=","))
dim(X)
library(readr)

Genenum <- read_csv("Data/Genenum.csv", col_names = FALSE)
colnames(X) <- unlist(Genenum)

# Standardise gene expression levels
X.std           <- scale(X)
# Remove predictor gene, which is tightly related to Ro131 as done in Segal et al
X.std           <- X.std[,c(1:6077,6079:6320)]

# Read in the response variable, entered manually
Y <- c(143,84,98,83,153,141,191,130,744,381,1047,806,621,849,
       475,966,708,487,1447,1693,1731,1025,376,126,102,149,91,
       153,235,68)
Y.std <- scale(Y)


n           <- dim(X.std)[1]
n # 30
p           <- dim(X.std)[2]
p #6319

type.order<- function(n){
  # max(2, floor((n)^(1/2)) - 3)
  max(1, floor(n^(1/3))-1)
}

set.seed(123)
mypval <- numeric()
# mypval_1 <- numeric()
pCorrs <- numeric()
dCorrs <- numeric()
for(i in 1:p){
  mypval[i] <- indeptest(X.std[,i],Y.std)$pvalue
  # mypval_1[i] <- indeptest(X[,i],Y, order = c(type.order(n),1))$pvalue
  pCorrs[i] <- abs(cor(X.std[,i],Y.std,method = "pearson"))
  dCorrs[i] <- abs(dcor(X.std[,i],y = Y.std))
  print(i)
}

pCorrOrder = order(pCorrs, decreasing = T)
dCorrOrder = order(dCorrs, decreasing = T)
BnCorrOrder = order(mypval)

d1 = floor(n/log(n))
d2 = 2*floor(n/log(n))
c(d1,d2)


data.frame(pCorrOrder[1:d1],
           dCorrOrder[1:d1],
           BnCorrOrder[1:d1])


mypval.BH <- p.adjust(mypval, method = "bonferroni")
length(which(mypval.BH<=0.05))
sel.index <- which(mypval.BH<=0.05)
colnames(X)[sel.index]

sel.index %in% BnCorrOrder[1:d1]
colnames(X)[BnCorrOrder[1:d1]]

############## Plot the scatterplot of the expression of z_Rol vs.
#############  the top d=d1 genes ranked by Bn-screen


require(graphics)
top=d1

pdf("/Volumes/GoogleDrive/Il mio Drive/autodep/R/All2/figures/cardio.pdf", width=8, height=12)

par(mfrow=c(4,2), mar=c(7,7,4,2)+.1)
for(i in 1:top){
  
  # scatter.smooth(XX[,CIS3ord.scaled[i]],YY, span = 0.6, degree = 1,
  #                   family = c("symmetric", "gaussian"),
  #                   xlab=paste("X [Rank ",i,"]"),
  #                   ylab="Y",
  #                   main=paste("CIS3.adj Rank", i),
  #                  ylim = range(YY,na.rm = TRUE),
  #                   evaluation = 100, lpars = list(col=2,lwd=2))
  Bn.loess=loess(Y.std ~ X.std[,BnCorrOrder[i]], 
                 span = 0.6, 
                 degree = 1, family = c("gaussian"),
                 evaluation = 100)
  Bn.loess.pred=predict(Bn.loess,se=T)
  ord=order(Bn.loess$x)
  
  
  plot(Bn.loess$x,Bn.loess$y,
       #xlab=paste(gene_name3[CIS3ord.scaled[i]],' (',probeIDsGenes3[CIS3ord.scaled[i],3],')'),
       xlab=colnames(X)[BnCorrOrder[i]],
       ylab="z_Rol", main=paste(c("Rank ", i),collapse = ""),
       font.lab=1.5,
       font.main=2,
       ylim = range(Y.std, na.rm = T),
       pch=21,
       cex.main=2,
       cex.lab=2,
       cex.axis=2,
       #xaxt="n",
       type="n")
  points(Bn.loess$x,Bn.loess$y, cex=2)
  lines(Bn.loess$x[ord],Bn.loess$fit[ord],col=2,lwd=2)
  lines(Bn.loess$x[ord],Bn.loess$fit[ord]+2*Bn.loess.pred$se.fit[ord],
        lty=2, lwd=2)
  lines(Bn.loess$x[ord],Bn.loess$fit[ord]-2*Bn.loess.pred$se.fit[ord],
        lty=2, lwd=2)
  
}

dev.off()

# same but using ggplot2
library(ggplot2)

dt <- data.frame()
for(i in 1:top){
  dt <- rbind(dt,
              data.frame(y = Y.std, x = X.std[, BnCorrOrder[i]],
                         gene = paste0(i, ". ", colnames(X.std)[BnCorrOrder[i]]))
              )  
}

ggplot(dt, aes(x=x, y=y)) +
  geom_point()+
  geom_smooth(method = "loess", formula = y~x, level = 0.95) +
  facet_wrap(vars(gene), 2, 4) +
  ylab("z_Rol") +
  xlab("gene") +
  theme_bw()

ggsave("cardio.pdf", width = 5, height = 3)

##### Zoom in the top d=d1 genes ranked by B-n screen ############## 

colnames(X)[BnCorrOrder[1:8]]

top.d1.3screens = list(
  PC=pCorrOrder[1:d1], 
  DC=dCorrOrder[1:d1],   
  Bn=BnCorrOrder[1:d1])




# mgcv --------------------------------------------------------------------
## no non viene bene : \
# 
# library(mgcv)
# 
# X.std.Bn.name <- colnames(X.std[,BnCorrOrder[1:d1]])
# 
# 
# Bn.mgcv <- gam(Y.std ~ 
#                  s(X.std[,BnCorrOrder[1]], bs = 'cr', k = 3) + s(X.std[,BnCorrOrder[2]], bs = 'cr', k = 3)+
#                  s(X.std[,BnCorrOrder[3]], bs = 'cr', k = 3) + s(X.std[,BnCorrOrder[4]], bs = 'cr', k = 3)+
#                  s(X.std[,BnCorrOrder[5]], bs = 'cr', k = 3) + s(X.std[,BnCorrOrder[6]], bs = 'cr', k = 3) +
#                  s(X.std[,BnCorrOrder[7]], bs = 'cr', k = 3) + s(X.std[,BnCorrOrder[8]], bs = 'cr', k = 3)
# )
# 
# summary(Bn.mgcv)$dev.expl
# summary(Bn.mgcv)$r.sq
# 
# PC.mgcv <- gam(Y.std ~ 
#                  s(X.std[,pCorrOrder[1]], bs = 'cr', k = 3) + s(X.std[,pCorrOrder[2]], bs = 'cr', k = 3)+
#                  s(X.std[,pCorrOrder[3]], bs = 'cr', k = 3) + s(X.std[,pCorrOrder[4]], bs = 'cr', k = 3)+
#                  s(X.std[,pCorrOrder[5]], bs = 'cr', k = 3) + s(X.std[,pCorrOrder[6]], bs = 'cr', k = 3) +
#                  s(X.std[,pCorrOrder[7]], bs = 'cr', k = 3) + s(X.std[,pCorrOrder[8]], bs = 'cr', k = 3)
# )
# 
# summary(PC.mgcv)$dev.expl
# summary(PC.mgcv)$r.sq
# 
# 
# DC.mgcv <- gam(Y.std ~ 
#                  s(X.std[,dCorrOrder[1]], bs = 'cr', k = 3) + s(X.std[,dCorrOrder[2]], bs = 'cr', k = 3)+
#                  s(X.std[,dCorrOrder[3]], bs = 'cr', k = 3) + s(X.std[,dCorrOrder[4]], bs = 'cr', k = 3)+
#                  s(X.std[,dCorrOrder[5]], bs = 'cr', k = 3) + s(X.std[,dCorrOrder[6]], bs = 'cr', k = 3) +
#                  s(X.std[,dCorrOrder[7]], bs = 'cr', k = 3) + s(X.std[,dCorrOrder[8]], bs = 'cr', k = 3)
# )
# 
# summary(DC.mgcv)$dev.expl
# summary(DC.mgcv)$r.sq


# Gamsel ------------------------------------------------------------------


set.seed(42)

nFolds = 5
foldid = sample(rep(seq(nFolds), 
                    length.out = nrow(X.std)))
            

#################################     Bn     #######################################

top=1500

X.std.Bn = X.std[, BnCorrOrder[1:top]]

Bn.cvGAMSEL.train =   cv.gamsel(X.std.Bn, Y.std,
                                family = "gaussian",
                                degrees = rep(5,top),
                                #gamma = gamsel.gamma,
                                # The default value is `gamma = 0.4`.
                                nfolds = nFolds,
                                foldid = foldid)

bestLambda.Bn = Bn.cvGAMSEL.train$lambda.1se
bestCVM.index = Bn.cvGAMSEL.train$index.1se
nzero.Train.Bn = Bn.cvGAMSEL.train$nzero[bestCVM.index]

Bn.gamsel.train = gamsel(X.std.Bn, Y.std,
                         family = "gaussian",
                         degrees = rep(5,top),
                         lambda = bestLambda.Bn)
devExpl.Bn = Bn.gamsel.train$dev.ratio
devExpl.Bn ## in-sample deviance explained 
# 0.9704428

nzero.Train.Bn 
# 28

#################################     PC     #######################################
X.std.PC = X.std[, pCorrOrder[1:top]]

PC.cvGAMSEL.train =   cv.gamsel(X.std.PC, Y.std,
                                family = "gaussian",
                                degrees = rep(5,top),
                                nfolds = nFolds,
                                foldid = foldid)

bestLambda.PC = PC.cvGAMSEL.train$lambda.1se
bestCVM.index = PC.cvGAMSEL.train$index.1se
nzero.Train.PC = PC.cvGAMSEL.train$nzero[bestCVM.index]

PC.gamsel.train = gamsel(X.std.PC, Y.std,
                         family = "gaussian",
                         degrees = rep(5,top),
                         lambda = bestLambda.PC)
devExpl.PC = PC.gamsel.train$dev.ratio
devExpl.PC
# 0.9490945
nzero.Train.PC
# 25


#################################     DC     #######################################
X.std.DC = X.std[, dCorrOrder[1:top]]

DC.cvGAMSEL.train =   cv.gamsel(X.std.DC, Y.std,
                                family = "gaussian",
                                degrees = rep(5,top),
                                nfolds = nFolds,
                                foldid = foldid)

bestLambda.DC = DC.cvGAMSEL.train$lambda.1se
bestCVM.index = DC.cvGAMSEL.train$index.1se
nzero.Train.DC = DC.cvGAMSEL.train$nzero[bestCVM.index]

DC.gamsel.train = gamsel(X.std.DC, Y.std,
                         family = "gaussian",
                         degrees = rep(5,top),
                         lambda = bestLambda.PC)
devExpl.DC = DC.gamsel.train$dev.ratio
devExpl.DC
# 0.9498063
nzero.Train.DC
# 29


#################################      GAMSEL      #######################################
## using all predictors
cvGAMSEL.train =   cv.gamsel(X.std, Y.std,
                             family = "gaussian",
                             degrees = rep(5,top),
                             nfolds = nFolds,
                             foldid = foldid)

bestLambda.gamsel = cvGAMSEL.train$lambda.1se
bestCVM.index = cvGAMSEL.train$index.1se
nzero.Train.GAMSEL = cvGAMSEL.train$nzero[bestCVM.index]

gamsel.train = gamsel(X.std, Y.std,
                      family = "gaussian",
                      degrees = rep(5,top),
                      lambda = bestLambda.gamsel)
devExpl.GAMSEL = gamsel.train$dev.ratio
devExpl.GAMSEL

# [1] 0.6952402
nzero.Train.GAMSEL
# 15

#################################      RAND     #######################################

nRep=200
devExpl.RAND = nzero.Train.RAND = rep(NA, nRep)


for(run in 1:nRep){
  
  print(run)
  
  RANDord = sample(1:ncol(X.std), top)
  
  XData.scaled.RAND = X.std[, RANDord]
  
  tryCatch({
    RAND.cvGAMSEL.train =   cv.gamsel(XData.scaled.RAND, Y.std,
                                      family = "gaussian",
                                      degrees = rep(5,top),
                                      nfolds = nFolds,
                                      foldid = foldid)
    
    bestLambda = RAND.cvGAMSEL.train$lambda.1se
    bestCVM.index = RAND.cvGAMSEL.train$index.1se
    nzero.Train.RAND[run] = RAND.cvGAMSEL.train$nzero[bestCVM.index]
    
    RAND.gamsel.train = gamsel(XData.scaled.RAND, Y.std,
                               family = "gaussian",
                               degrees = rep(5,top),
                               lambda = bestLambda)
    devExpl.RAND[run] = RAND.gamsel.train$dev.ratio
  }, error=function(e){})
}


boxplot(devExpl.RAND)

summary(devExpl.RAND)
mean(devExpl.RAND)
# 0.5220799
median(devExpl.RAND)
# 0.536852

my.summary<-data.frame(dv.exp =c(
  devExpl.Bn*100,
  devExpl.PC*100,
  devExpl.DC*100,
  devExpl.GAMSEL*100, median(devExpl.RAND)*100), 
  n.nzeros =c(
    nzero.Train.Bn,
    nzero.Train.PC,
    nzero.Train.DC,
    nzero.Train.GAMSEL, 
    median(nzero.Train.RAND)
  ))


my.summary <- t(my.summary)

colnames(my.summary) <- c("Bn", "PC", "DC", "GAMSEL", "RAND")

library(xtable)
xtable(my.summary)

# intersection ------------------------------------------------------------
top.8.3screens = list(Bn=BnCorrOrder[1:d1], 
                      PC=pCorrOrder[1:d1],
                      DC=dCorrOrder[1:d1])

top.1500.3screens  = list(Bn=BnCorrOrder[1:top], 
                          PC=pCorrOrder[1:top],
                          DC=dCorrOrder[1:top])

common.3screens = matrix(NA, 3, 3)
colnames(common.3screens) = rownames(common.3screens) = 
  c("Bn", "PC", "DC")

common.3screens_prop = matrix(NA, 3, 3)
colnames(common.3screens_prop) = rownames(common.3screens_prop) = 
  c("Bn", "PC", "DC")

for(i in 1:2){
  for(j in (i+1):3){
    
    common.3screens[i,j] = length(intersect(unlist(top.8.3screens[i]), unlist(top.8.3screens[j])))
    common.3screens[j,i] = length(intersect(unlist(top.1500.3screens[i]), unlist(top.1500.3screens[j])))
    
    common.3screens_prop[i,j] = length(intersect(unlist(top.8.3screens[i]), unlist(top.8.3screens[j])))/d1
    common.3screens_prop[j,i] = length(intersect(unlist(top.1500.3screens[i]), unlist(top.1500.3screens[j])))/top
  }
}

common.3screens
round(common.3screens_prop,3)


xtable(common.3screens)
xtable(common.3screens_prop)


Reduce(intersect, list(BnCorrOrder[1:d1],
                       pCorrOrder[1:d1],
                       dCorrOrder[1:d1])) #3655  349
which(BnCorrOrder[1:d1] %in% c(3655, 349)) # 1, 2

colnames(X)[c(3655, 349)]

Reduce(intersect, list(BnCorrOrder[1:d1],
                       pCorrOrder[1:d1])) #3655  349


Reduce(intersect, list(BnCorrOrder[1:d1],
                       dCorrOrder[1:d1])) #3655  349 1107
which(BnCorrOrder[1:d1] %in% c(3655,  349, 1107))  # 1, 2, 6

colnames(X.std)[1107]

# out of sample prediction error ------------------------------------------

library(tictoc)


nRep=200
train.Prop = 0.7


## Bn -----------------------------------------------------

devExpl.Bn.run=rmse.Bn.run=nzero.Train.Bn.run=rep(NA, nRep)
nTrainIDs=matrix(NA,nRep,ceiling(nrow(X.std)*train.Prop))

nFolds=5

library(tictoc)
tic()

for(run in 1:nRep){

  print(run)
  set.seed(run)

  nTrain = sample(1:nrow(X.std), ceiling(nrow(X.std)*train.Prop))
  nTrainIDs[run,] = nTrain

  foldid = sample(rep(seq(nFolds), length.out = length(nTrain)))

  yData_Train = Y.std[nTrain]
  yData_Test  = Y.std[-nTrain]

  XData_Train_Bn = X.std[nTrain, BnCorrOrder[1:top]]
  XData_Test_Bn  = X.std[-nTrain, BnCorrOrder[1:top]]

  Bn.cvGAMSEL.train.run =   cv.gamsel(XData_Train_Bn, yData_Train,
                                    family = "gaussian",
                                    degrees = rep(5,top),
                                    nfolds = nFolds,
                                    foldid = foldid)

  bestLambda.run = Bn.cvGAMSEL.train.run$lambda.1se
  bestCVM.index.run = Bn.cvGAMSEL.train.run$index.1se
  nzero.Train.Bn.run[run] = Bn.cvGAMSEL.train.run$nzero[bestCVM.index.run]

  Bn.gamsel.train.run = gamsel(XData_Train_Bn, yData_Train,
                             family = "gaussian",
                             degrees = rep(5,top),
                             lambda = bestLambda.run)
  devExpl.Bn.run[run] = Bn.gamsel.train.run$dev.ratio

  Bn.gamsel.test = predict(Bn.gamsel.train.run, XData_Test_Bn,
                                    type = "response")
  rmse.Bn.run[run] = sqrt(mean((yData_Test - Bn.gamsel.test)^2))
}

toc() 
# 1281.377 sec elapsed


mylist.Bn.run <- list(devExpl.Bn.run = devExpl.Bn.run, 
                    rmse.Bn.run = rmse.Bn.run, 
                    nzero.Train.Bn.run  = nzero.Train.Bn.run)
save(mylist.Bn.run, file="results/Bn_run_nRep200.RData")


# KCV estimate 
mean(rmse.Bn.run)
# 0.86038
boxplot(rmse.Bn.run)

mean(devExpl.Bn.run)

## PC -----------------------------------------------------

devExpl.PC.run=rmse.PC.run=nzero.Train.PC.run=rep(NA, nRep)
nTrainIDs=matrix(NA,nRep,ceiling(nrow(X.std)*train.Prop))

nFolds=5

tic()
for(run in 1:nRep){
  
  print(run)
  set.seed(run)
  
  nTrain = sample(1:nrow(X.std), ceiling(nrow(X.std)*train.Prop))
  nTrainIDs[run,] = nTrain
  
  foldid = sample(rep(seq(nFolds), length.out = length(nTrain)))
  
  yData_Train = Y.std[nTrain]
  yData_Test  = Y.std[-nTrain]
  
  XData_Train_PC = X.std[nTrain, pCorrOrder[1:top]]
  XData_Test_PC  = X.std[-nTrain, pCorrOrder[1:top]]
  
  PC.cvGAMSEL.train.run =   cv.gamsel(XData_Train_PC, yData_Train,
                                      family = "gaussian",
                                      degrees = rep(5,top),
                                      nfolds = nFolds,
                                      foldid = foldid)
  
  bestLambda.run = PC.cvGAMSEL.train.run$lambda.1se
  bestCVM.index.run = PC.cvGAMSEL.train.run$index.1se
  nzero.Train.PC.run[run] = PC.cvGAMSEL.train.run$nzero[bestCVM.index.run]
  
  PC.gamsel.train.run = gamsel(XData_Train_PC, yData_Train,
                               family = "gaussian",
                               degrees = rep(5,top),
                               lambda = bestLambda.run)
  devExpl.PC.run[run] = PC.gamsel.train.run$dev.ratio
  
  PC.gamsel.test = predict(PC.gamsel.train.run, XData_Test_PC,
                           type = "response")
  rmse.PC.run[run] = sqrt(mean((yData_Test - PC.gamsel.test)^2))
}
toc() 


mylist.PC.run <- list(devExpl.PC.run = devExpl.PC.run, 
                      rmse.PC.run = rmse.PC.run, 
                      nzero.Train.PC.run  = nzero.Train.PC.run)
save(mylist.PC.run, file="results/PC_run_nRep200.RData")



# KCV estimate 
mean(rmse.PC.run)
# [1] 0.7799849

boxplot(rmse.Bn.run, rmse.PC.run)

mean(devExpl.PC.run)

## DC -----------------------------------------------------

devExpl.DC.run=rmse.DC.run=nzero.Train.DC.run=rep(NA, nRep)
nTrainIDs=matrix(NA,nRep,ceiling(nrow(X.std)*train.Prop))

nFolds=5

tic()
for(run in 1:nRep){
  
  print(run)
  set.seed(run)
  
  nTrain = sample(1:nrow(X.std), ceiling(nrow(X.std)*train.Prop))
  nTrainIDs[run,] = nTrain
  
  foldid = sample(rep(seq(nFolds), length.out = length(nTrain)))
  
  yData_Train = Y.std[nTrain]
  yData_Test  = Y.std[-nTrain]
  
  XData_Train_DC = X.std[nTrain, dCorrOrder[1:top]]
  XData_Test_DC  = X.std[-nTrain, dCorrOrder[1:top]]
  
  DC.cvGAMSEL.train.run =   cv.gamsel(XData_Train_DC, yData_Train,
                                      family = "gaussian",
                                      degrees = rep(5,top),
                                      nfolds = nFolds,
                                      foldid = foldid)
  
  bestLambda.run = DC.cvGAMSEL.train.run$lambda.1se
  bestCVM.index.run = DC.cvGAMSEL.train.run$index.1se
  nzero.Train.DC.run[run] = DC.cvGAMSEL.train.run$nzero[bestCVM.index.run]
  
  DC.gamsel.train.run = gamsel(XData_Train_DC, yData_Train,
                               family = "gaussian",
                               degrees = rep(5,top),
                               lambda = bestLambda.run)
  devExpl.DC.run[run] = DC.gamsel.train.run$dev.ratio
  
  DC.gamsel.test = predict(DC.gamsel.train.run, XData_Test_DC,
                           type = "response")
  rmse.DC.run[run] = sqrt(mean((yData_Test - DC.gamsel.test)^2))
}
toc() 

mylist.DC.run <- list(devExpl.DC.run = devExpl.DC.run, 
                      rmse.DC.run = rmse.DC.run, 
                      nzero.Train.DC.run  = nzero.Train.DC.run)
save(mylist.DC.run, file="results/DC_run_nRep200.RData")


# KCV estimate 
mean(rmse.DC.run)
# 0.8191413

## GAMSEL -----------------------------------------------------

devExpl.GAMSEL.run=rmse.GAMSEL.run=nzero.Train.GAMSEL.run=rep(NA, nRep)
nTrainIDs=matrix(NA,nRep,ceiling(nrow(X.std)*train.Prop))

nFolds=5

tic()
for(run in 1:nRep){
  
  print(run)
  set.seed(run)
  
  nTrain = sample(1:nrow(X.std), ceiling(nrow(X.std)*train.Prop))
  nTrainIDs[run,] = nTrain
  
  foldid = sample(rep(seq(nFolds), length.out = length(nTrain)))
  
  yData_Train = Y.std[nTrain]
  yData_Test  = Y.std[-nTrain]
  
  XData_Train_GAMSEL = X.std[nTrain, ]
  XData_Test_GAMSEL  = X.std[-nTrain, ]
  
  GAMSEL.cvGAMSEL.train.run =   cv.gamsel(XData_Train_GAMSEL, yData_Train,
                                      family = "gaussian",
                                      degrees = rep(5,top),
                                      nfolds = nFolds,
                                      foldid = foldid)
  
  bestLambda.run = GAMSEL.cvGAMSEL.train.run$lambda.1se
  bestCVM.index.run = GAMSEL.cvGAMSEL.train.run$index.1se
  nzero.Train.GAMSEL.run[run] = GAMSEL.cvGAMSEL.train.run$nzero[bestCVM.index.run]
  
  GAMSEL.gamsel.train.run = gamsel(XData_Train_GAMSEL, yData_Train,
                               family = "gaussian",
                               degrees = rep(5,top),
                               lambda = bestLambda.run)
  devExpl.GAMSEL.run[run] = GAMSEL.gamsel.train.run$dev.ratio
  
  GAMSEL.gamsel.test = predict(GAMSEL.gamsel.train.run, XData_Test_GAMSEL,
                           type = "response")
  rmse.GAMSEL.run[run] = sqrt(mean((yData_Test - GAMSEL.gamsel.test)^2))
}
toc() 


# KCV estimate 
mean(rmse.GAMSEL.run)
#

mylist.GAMSEL.run <- list(devExpl.GAMSEL.run = devExpl.GAMSEL.run, 
                      rmse.GAMSEL.run = rmse.GAMSEL.run, 
                      nzero.Train.GAMSEL.run  = nzero.Train.GAMSEL.run)
save(mylist.GAMSEL.run, file="results/GAMSEL_run_nRep200.RData")



## RANDSEL -----------------------------------------------------


devExpl.RAND=rmse.RAND=nzero.Train.RAND=rep(NA, nRep)
nTrainIDs=matrix(NA,nRep,ceiling(nrow(X.std)*train.Prop))


tic()
for(run in 1:nRep){
  
  print(run)
  set.seed(run)
  
  RANDord = sample(1:ncol(X.std), top)
  
  nTrain = sample(1:nrow(X.std), ceiling(nrow(X.std)*train.Prop))
  nTrainIDs[run,] = nTrain
  
  foldid = sample(rep(seq(nFolds), length.out = length(nTrain)))
  
  yData_Train = Y.std[nTrain]
  yData_Test = Y.std[-nTrain]
  
  XData_Train_RAND = X.std[nTrain, RANDord]
  XData_Test_RAND = X.std[-nTrain, RANDord]
  
  RAND.cvGAMSEL.train =   cv.gamsel(XData_Train_RAND, yData_Train,
                                    family = "gaussian",
                                    degrees = rep(5,top),
                                    nfolds = nFolds,
                                    foldid = foldid)
  
  bestLambda = RAND.cvGAMSEL.train$lambda.1se
  bestCVM.index = RAND.cvGAMSEL.train$index.1se
  nzero.Train.RAND[run] = RAND.cvGAMSEL.train$nzero[bestCVM.index]
  
  RAND.gamsel.train = gamsel(XData_Train_RAND, yData_Train,
                             family = "gaussian",
                             degrees = rep(5,top),
                             lambda = bestLambda)
  devExpl.RAND[run] = RAND.gamsel.train$dev.ratio
  
  RAND.gamsel.test = predict(RAND.gamsel.train, XData_Test_RAND,
                                    type = "response")
  RAND.gamsel.yTest = RAND.gamsel.test
  rmse.RAND[run] = sqrt(mean((yData_Test - RAND.gamsel.yTest)^2))
}
toc()

mylist.rand <- list(devExpl.RAND = devExpl.RAND, 
       rmse.RAND = rmse.RAND, 
       nzero.Train.RAND  = nzero.Train.RAND)
save(mylist.rand, file="results/RANDgamsel_nRep200.RData")



# reading lists  ----------------------------------------------------------

load("results/Bn_run_nRep200.RData")
load("results/PC_run_nRep200.RData")
load("results/DC_run_nRep200.RData")
load("results/GAMSEL_run_nRep200.RData")
load("results/RANDgamsel_nRep200.RData")

pdf("figures/cardio_rmse.pdf", width=8, height=5)
boxplot(mylist.Bn.run$rmse.Bn.run, 
        mylist.PC.run$rmse.PC.run, 
        mylist.DC.run$rmse.DC.run,  
        mylist.GAMSEL.run$rmse.GAMSEL.run,
        mylist.rand$rmse.RAND, 
        ylim=c(0,4), col="white", xlab="Methods",
        ylab="RMSPE (30% validation sets)", 
        names=c("Bn", "PC", "DC", "GAMSEL", "RAND"))
dev.off()

pdf("figures/cardio_devExpl.pdf", width=8, height=5)
boxplot(mylist.Bn.run$devExpl.Bn.run * 100, 
        mylist.PC.run$devExpl.PC.run* 100, 
        mylist.DC.run$devExpl.DC.run* 100, 
        mylist.GAMSEL.run$devExpl.GAMSEL* 100, 
        mylist.rand$devExpl.RAND* 100, col="white", xlab="Methods",
        ylab="Percent deviance explained (70% training sets)", 
        names=c("Bn", "PC", "DC", "GAMSEL", "RAND"))
dev.off()


pdf("figures/cardio_nonzeros.pdf", width=8, height=5)
boxplot(mylist.Bn.run$nzero.Train.Bn.run, 
        mylist.PC.run$nzero.Train.PC.run, 
        mylist.DC.run$nzero.Train.DC.run, 
        mylist.GAMSEL.run$nzero.Train.GAMSEL.run, 
        mylist.rand$nzero.Train.RAND, col="white", xlab="Methods",
        ylab="Number of non-zero coefficient estimates (70% training sets)", 
        names=c("Bn", "PC", "DC", "GAMSEL", "RAND"))
dev.off()

c(mean(devExpl.Bn.run), mean(devExpl.PC.run),
  mean(devExpl.DC.run), mean(devExpl.GAMSEL), 
  mean(devExpl.RAND))


# the same but using ggplot2
df_exdev <- data.frame(value = c(mylist.Bn.run$devExpl.Bn.run, 
                                 mylist.PC.run$devExpl.PC.run, 
                                 mylist.DC.run$devExpl.DC.run,
                                 mylist.GAMSEL.run$devExpl.GAMSEL.run,
                                 mylist.rand$devExpl.RAND)*100,
                       method = rep(c("Bn", "PC", "DC", "GAMSEL", "RANDOM"), each = 200),
                       measure = "Explained deviance on train set (%)"
                       )

df_rmspe <- data.frame(value = c(mylist.Bn.run$rmse.Bn.run, 
                                 mylist.PC.run$rmse.PC.run, 
                                 mylist.DC.run$rmse.DC.run,
                                 mylist.GAMSEL.run$rmse.GAMSEL.run,
                                 mylist.rand$rmse.RAND),
                       method = rep(c("Bn", "PC", "DC", "GAMSEL", "RANDOM"), each = 200),
                       measure = "RMSPE on validation set"
)

df_nzero <- data.frame(value = c(mylist.Bn.run$nzero.Train.Bn.run, 
                                 mylist.PC.run$nzero.Train.PC.run, 
                                 mylist.DC.run$nzero.Train.DC.run,
                                 mylist.GAMSEL.run$nzero.Train.GAMSEL.run,
                                 mylist.rand$nzero.Train.RAND),
                       method = rep(c("Bn", "PC", "DC", "GAMSEL", "RANDOM"), each = 200),
                       measure = "Number of non-zero coefficients (%)"
)

df <- rbind(df_exdev, df_nzero, df_rmspe)

noplot <- (df$measure == "RMSPE on validation set") & (df$value > 2)

df$method <- factor(df$method, levels = c("Bn", "PC", "DC", "GAMSEL", "RANDOM"))

ggplot(df[!noplot, ], aes(x = method, y = value)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap(vars(measure), nrow = 3, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme_bw()

ggsave("cardio_boxplots.pdf", width = 5, height = 5)

# LOOCV -------------------------------------------------------------------

devExpl.Bn=Bn_oneout=nzero.Train.Bn=rep(NA, n)

devExpl.PC=PC_oneout=nzero.Train.PC=rep(NA, n)

devExpl.DC=DC_oneout=nzero.Train.DC=rep(NA, n)

set.seed(42)

for (i in 1:n){
  
  Y.std_Train = Y.std[-i]
  Y.std_Test  = Y.std[i]
  
  foldid = sample(rep(seq(nFolds), length.out = n-1))
  
  ## Bn --------------------------------------------------------
  
  Xstd_Train_Bn = X.std[-i, BnCorrOrder[1:top]]
  Xstd_Test_Bn  = X.std[i, BnCorrOrder[1:top]]
  
  Bn.cvGAMSEL.train =   cv.gamsel(Xstd_Train_Bn, Y.std_Train,
                                  family = "gaussian",
                                  degrees = rep(5,top),
                                  #gamma = gamsel.gamma,
                                  nfolds = nFolds,
                                  foldid = foldid)
  bestLambda.Bn = Bn.cvGAMSEL.train$lambda.1se
  bestCVM.index = Bn.cvGAMSEL.train$index.1se
  nzero.Train.Bn[i] = Bn.cvGAMSEL.train$nzero[bestCVM.index]
  
  Bn.gamsel.train.run = gamsel(Xstd_Train_Bn, Y.std_Train,
                           family = "gaussian",
                           degrees = rep(5,top),
                           #gamma = gamsel.gamma,
                           lambda = bestLambda.Bn)
  devExpl.Bn[i] = Bn.gamsel.train.run$dev.ratio
  
  yhat_i = predict(Bn.gamsel.train.run, 
                           t(as.matrix(Xstd_Test_Bn)),
                           type = "response")
  Bn_oneout[i] <- ( Y.std_Test -  yhat_i )^2
  
  
  ## PC --------------------------------------------------------
  
  Xstd_Train_PC = X.std[-i, pCorrOrder[1:top]]
  Xstd_Test_PC  = X.std[i, pCorrOrder[1:top]]
  
  PC.cvGAMSEL.train =   cv.gamsel(Xstd_Train_PC, Y.std_Train,
                                  family = "gaussian",
                                  degrees = rep(5,top),
                                  #gamma = gamsel.gamma,
                                  nfolds = nFolds,
                                  foldid = foldid)
  bestLambda.PC = PC.cvGAMSEL.train$lambda.1se
  bestCVM.index = PC.cvGAMSEL.train$index.1se
  nzero.Train.PC[i] = PC.cvGAMSEL.train$nzero[bestCVM.index]
  
  PC.gamsel.train.run = gamsel(Xstd_Train_PC, Y.std_Train,
                               family = "gaussian",
                               degrees = rep(5,top),
                               #gamma = gamsel.gamma,
                               lambda = bestLambda.PC)
  devExpl.PC[i] = PC.gamsel.train.run$dev.ratio
  
  yhat_i = predict(PC.gamsel.train.run, 
                   t(as.matrix(Xstd_Test_PC)),
                   type = "response")
  PC_oneout[i] <- ( Y.std_Test -  yhat_i )^2
  
  
  ## DC --------------------------------------------------------
  
  Xstd_Train_DC = X.std[-i, dCorrOrder[1:top]]
  Xstd_Test_DC  = X.std[i, dCorrOrder[1:top]]
  
  DC.cvGAMSEL.train =   cv.gamsel(Xstd_Train_DC, Y.std_Train,
                                  family = "gaussian",
                                  degrees = rep(5,top),
                                  #gamma = gamsel.gamma,
                                  nfolds = nFolds,
                                  foldid = foldid)
  bestLambda.DC = DC.cvGAMSEL.train$lambda.1se
  bestCVM.index = DC.cvGAMSEL.train$index.1se
  nzero.Train.DC[i] = DC.cvGAMSEL.train$nzero[bestCVM.index]
  
  DC.gamsel.train.run = gamsel(Xstd_Train_DC, Y.std_Train,
                               family = "gaussian",
                               degrees = rep(5,top),
                               #gamma = gamsel.gamma,
                               lambda = bestLambda.DC)
  devExpl.DC[i] = DC.gamsel.train.run$dev.ratio
  
  yhat_i = predict(DC.gamsel.train.run, 
                   t(as.matrix(Xstd_Test_DC)),
                   type = "response")
  DC_oneout[i] <- ( Y.std_Test -  yhat_i )^2
  
  print(i)
}

mean(Bn_oneout)
# [1] 0.3621615

mean(PC_oneout)
#[1] 0.3478313

mean(DC_oneout)
# [1] 0.3740485

boxplot(devExpl.Bn, devExpl.PC, devExpl.DC)
# 0.6845866

mean(devExpl.PC)
# 0.7233857

mean(devExpl.DC)
# 0.7220343
