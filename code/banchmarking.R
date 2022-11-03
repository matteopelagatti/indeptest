library(dplyr)
library(ggplot2)
source("indeptest_with_cat.R")
source("libraries.R")
source("utilities.R")


# last august 2022 --------------------------------------------------------

set.seed(123)

ss <- c(100, 200, 400, 800, 1600, 3200, 6400, 12800)           # sample size
ban1 <- vector("list", length(ss))    # will contain microbanchmark results

for(i in 1:length(ss)){
  cat("Simulating n =", ss[i], "\n")
  x <- rnorm(ss[i])
  y <- rnorm(ss[i])
  
  ban1[[i]] <- microbenchmark(
    "Bn"    = indeptest(x, y),
    "Hoef"  = hoeffd(x, y),
    "MIC"   = mine_stat(x, y),
    "dCov"  = dcov.test(x, y),
    "HHG"   = hhg(x, y),
    times = 100
  )
}

names(ban1) <- paste0("n = ", ss)

dt <- Reduce(rbind, ban1)
dt$n <- rep(ss, each = 500)

write.csv(dt, "exec_times.csv", row.names = FALSE)
dt <- read.csv("exec_times.csv")

dtsum <- dt |> group_by(expr, n) |> summarise(median = median(time),
                                              p05 = quantile(time, 0.05),
                                              p95 = quantile(time, 0.95))

dtsum <- dtsum |> rename(test = expr)

dtsum |> ggplot(aes(x = n, y = median/1000000000,
                    ymin = p05/1000000000,
                    ymax = p95/1000000000,
                    color = test, lty = test, shape = test, fill = test)) +
  geom_line() +
  # geom_ribbon(alpha = 0.1)+
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  ylab("seconds") +
  theme_bw()
ggsave("exec_times.pdf", width = 7, height = 3)



dtsum |> mutate(sec = median/1000000000) |> select(test, n, sec) |> View()
