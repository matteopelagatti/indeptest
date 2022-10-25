library(Hmisc)
library(mvtnorm)
library(ggplot2)
source("indeptest_with_cat.R")

nsim <- 10000
n <- c(25, 50, 100, 200)
r <- seq(0, .9, by = .1)


cont.table <- function(n, order, rho=0.5){
  mu <- c(0,0)
  Sigma <- matrix(c(1,rho,rho,1), ncol=2)
  x <- rmvnorm(n, mean=mu, sigma=Sigma)
  if(order == 3){
    x.ord <- as.factor(ifelse(x[,1]<= -.6, 1, ifelse(x[,1]>.6, 3, 2)))
    y.ord <- as.factor(ifelse(x[,2]<= -.6, 1, ifelse(x[,2]>.6, 3, 2)))
  }
  else if(order == 4){
    x.ord <- as.factor(ifelse(x[,1]<= -.8, 1,
                              ifelse(x[,1]>.8, 4, 
                                     ifelse(x[,1]>0 & x[,1]<=.8, 3, 2))))
    y.ord <- as.factor(ifelse(x[,2]<= -.8, 1,
                              ifelse(x[,2]>.8, 4, 
                                     ifelse(x[,2]>0 & x[,2]<=.8, 3, 2))))
  }
  else if(order == 5){
    x.ord <- as.factor(ifelse(x[,1]<= -1, 1,
                              ifelse(x[,1]>1, 5, 
                                     ifelse(x[,1]>-1 & x[,1]<=-.3, 2, 
                                            ifelse(x[1,]>-.3 & x[1,]<=0,3,4)))))
    y.ord <- as.factor(ifelse(x[,2]<= -1, 1,
                              ifelse(x[,2]>1, 5, 
                                     ifelse(x[,2]>-1 & x[,2]<=-.3, 2, 
                                            ifelse(x[2,]>-.3 & x[2,]<=0,3,4)))))
  }
  rbind(x.ord, y.ord)
}


set.seed(42)
mytest.pval3 <- array(NA, c(nsim, length(n), length(r)))
chisqtest.pval3 <- array(NA, c(nsim, length(n), length(r)))
hoeffd.pval3 <- array(NA, c(nsim, length(n), length(r)))
mytest.pval5 <- array(NA, c(nsim, length(n), length(r)))
chisqtest.pval5 <- array(NA, c(nsim, length(n), length(r)))
hoeffd.pval5 <- array(NA, c(nsim, length(n), length(r)))

for(k in 1: length(n)){
  for(j in 1:length(r)){
    for (i in 1:nsim){
      xy <- cont.table(n[k], order = 3, rho= r[j])
      mytest.pval3[i,k,j] <- indeptest(as.factor(xy[1, ]), as.factor(xy[2, ]), 
                                       basis ="dummy")$pvalue
      chisqtest.pval3[i,k,j] <- chisq.test(as.factor(xy[1, ]), as.factor(xy[2, ]))$p.value
      hoeffd.pval3[i,k,j] <- hoeffd(as.factor(xy[1, ]), as.factor(xy[2, ]))$P[1,2]
      xy <- cont.table(n[k], order = 5, rho= r[j])
      mytest.pval5[i,k,j] <- indeptest(as.factor(xy[1, ]), as.factor(xy[2, ]), 
                                       basis ="dummy")$pvalue
      chisqtest.pval5[i,k,j] <- chisq.test(as.factor(xy[1, ]), as.factor(xy[2, ]))$p.value
      hoeffd.pval5[i,k,j] <- hoeffd(as.factor(xy[1, ]), as.factor(xy[2, ]))$P[1,2]
    }
  }
  cat("i=",i, "j=", j, "k=", k, "\n")
}

dt <- data.frame(
  rejection_rate = as.numeric(colMeans(mytest.pval3 < 0.05, na.rm = TRUE)),
  n = n,
  rho = rep(r, each = length(n)),
  table = "3 x 3",
  test = "Bn"
)

dt <- rbind(dt,
            data.frame(
              rejection_rate = as.numeric(colMeans(chisqtest.pval3 < 0.05, na.rm = TRUE)),
              n = n,
              rho = rep(r, each = length(n)),
              table = "3 x 3",
              test = "ChiSqr"
            )
)

dt <- rbind(dt,
            data.frame(
              rejection_rate = as.numeric(colMeans(hoeffd.pval3 < 0.05, na.rm = TRUE)),
              n = n,
              rho = rep(r, each = length(n)),
              table = "3 x 3",
              test = "Hoef"
            )
)

dt <- rbind(dt,
            data.frame(
              rejection_rate = as.numeric(colMeans(mytest.pval5 < 0.05, na.rm = TRUE)),
              n = n,
              rho = rep(r, each = length(n)),
              table = "5 x 5",
              test = "Bn"
            )
)

dt <- rbind(dt,
            data.frame(
              rejection_rate = as.numeric(colMeans(chisqtest.pval5 < 0.05, na.rm = TRUE)),
              n = n,
              rho = rep(r, each = length(n)),
              table = "5 x 5",
              test = "ChiSqr"
            )
)

dt <- rbind(dt,
            data.frame(
              rejection_rate = as.numeric(colMeans(hoeffd.pval5 < 0.05, na.rm = TRUE)),
              n = n,
              rho = rep(r, each = length(n)),
              table = "5 x 5",
              test = "Hoef"
            )
)

data.table::fwrite(dt, "rejection_rates_nominal.csv")

dt$n <- paste0("n = ", dt$n)
dt$n <- factor(dt$n, levels = c("n = 25", "n = 50", "n = 100", "n = 200"))

dt %>% ggplot(aes(x = rho, y = rejection_rate, color = test, linetype = test, shape = test)) +
  facet_grid(vars(n), vars(table)) +
  geom_line() + geom_point() + geom_hline(yintercept = 0.05, lty = 2) +
  ylab("rejection rate") + xlab(expression(rho)) +
  theme_bw()
ggsave("power_nominal_plot.pdf", width = 7,height = 7)

