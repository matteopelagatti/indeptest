library(data.table)
library(HHG)
library(energy)
library(minerva)
library(Hmisc)
source("indeptest_with_cat.R")
source("noise.distributions.R")

# setup
n <- 100
nsim <- 10000

types <- 1:10
noises <- seq(0.1, 2, 0.1)

type_names <- c("Linear", "Quadratic", "Cubic", "Sine",
                "Fourth root", "Circle", "Two curves", "X", "Diamond", "XsdY")

# ---- generate and save simulated observations ----
# for every simulation the columns x, y and a permutation of y are stored
# the last column is used to compute the empirical size of the test
set.seed(202210)
for (ty in types) {
  for (no in noises) {
    sims <- matrix(nrow = n*nsim, ncol = 3, dimnames = list(NULL, c("x", "y", "yperm")))
    sims[, 1:2] <- t(datagen.noise(n*nsim, ty, no))
    sims[, 3] <- sims[sample.int(nrow(sims)), 2]
    fname <- paste0("sim_type", ty, "_noise", no, ".csv")
    fwrite(sims, fname)
  }
}

# ---- tests ----
# we compute the test statistics on (x, y) and store in *_dep variables
# we compute the test statistics on (x, permutated) and store in *_ind variables
# we store also the execution times
hhg_dep <- numeric(nsim)
hhg_ind <- numeric(nsim)
hhg_tim <- numeric(nsim)
dcov_dep <- numeric(nsim)
dcov_ind <- numeric(nsim)
dcov_tim <- numeric(nsim)
mic_dep <- numeric(nsim)
mic_ind <- numeric(nsim)
mic_tim <- numeric(nsim)
hoef_dep <- numeric(nsim)
hoef_ind <- numeric(nsim)
hoef_tim <- numeric(nsim)
our_dep <- numeric(nsim)
our_ind <- numeric(nsim)
our_tim <- numeric(nsim)
for (ty in types) {
  for (no in noises) {
    fname <- paste0("sim_type", ty, "_noise", no, ".csv")
    sims  <- as.matrix(fread(fname))
    cat("Working on", fname, "\n")
    for (i in 1:nsim) {
      X <- sims[((i-1)*n+1):(i*n), ]
      # HHG
      start_time <- Sys.time()
        Dx = as.matrix(dist((X[, 1]), diag = TRUE, upper = TRUE))
        Dy = as.matrix(dist((X[, 2]), diag = TRUE, upper = TRUE))
        Dy0 = as.matrix(dist((X[, 3]), diag = TRUE, upper = TRUE))
        hhg_dep[i] <- hhg.test(Dx, Dy,  nr.perm = 0)[[1]]
        hhg_ind[i] <- hhg.test(Dx, Dy0, nr.perm = 0)[[1]]
      end_time <- Sys.time()
      hhg_tim[i] <- end_time - start_time
      # dCov
      start_time <- Sys.time()
        dcov_dep[i] <- dcov.test(X[, 1], X[, 2])$statistic
        dcov_ind[i] <- dcov.test(X[, 1], X[, 3])$statistic
      end_time <- Sys.time()
      dcov_tim[i] <- end_time - start_time
      # MIC
      start_time <- Sys.time()
        mic_dep[i] <- mine_stat(X[, 1], X[, 2])
        mic_ind[i] <- mine_stat(X[, 1], X[, 3])
      end_time <- Sys.time()
      mic_tim[i] <- end_time - start_time
      # Hoeffding
      start_time <- Sys.time()
        hoef_dep[i] <- hoeffd(X[, 1], X[, 2])$D[2]
        hoef_ind[i] <- hoeffd(X[, 1], X[, 3])$D[2]
      end_time <- Sys.time()
      hoef_tim[i] <- end_time - start_time
      # Our test
      start_time <- Sys.time()
        our_dep[i] <- indeptest(X[, 1], X[, 2],
                                order = if (ty <= 5) c(3, 1) else c(3, 3))$stat
        our_ind[i] <- indeptest(X[, 1], X[, 3],
                                order = if (ty <= 5) c(3, 1) else c(3, 3))$stat
      end_time <- Sys.time()
      our_tim[i] <- end_time - start_time
    }
    tests <- data.frame(
      power = c(
        mean(hhg_dep  > quantile(hhg_ind,  0.95)),
        mean(dcov_dep > quantile(dcov_ind, 0.95)),
        mean(mic_dep  > quantile(mic_ind,  0.95)),
        mean(hoef_dep > quantile(hoef_ind, 0.95)),
        mean(our_dep  > quantile(our_ind,  0.95))
      ),
      time = sapply(list(hhg_tim, dcov_tim, mic_tim, hoef_tim, our_tim), mean),
      noise = no,
      type = type_names[ty],
      test = c("HHG", "dCov", "MIC", "Hoef", "Bn")
    )
    if (ty == types[1] & no == noises[1]) {
      fwrite(tests, "test_power.csv")
    } else {
      fwrite(tests, "test_power.csv", append = TRUE)
    }
  }
}

# plotting
library(ggplot2)
dt <- fread("test_power.csv")
dt$type <- factor(dt$type, levels = type_names)
dt |> ggplot(aes(x = noise, y = power, shape = test, linetype = test, color = test)) +
  geom_line() +
  geom_point() +
  facet_wrap(~type, ncol = 2) +
  theme_bw()
ggsave("power_plot.pdf", width = 7, height = 9)


# ranking the tests
library(dplyr)
# mean power
dt |> group_by(test) |> summarise(mean_power = mean(power)) |> arrange(desc(mean_power))
# mean rank
wdt <- dt |> select(-time) |> tidyr::pivot_wider(names_from = test, values_from = power)
rdt <- t(apply(as.matrix(wdt[, 3:7]), 1, frankv, order = -1L, ties.method = "min"))
colnames(rdt) <- paste0("rank_", names(wdt[, 3:7]))
rankdt <- cbind(wdt, rdt)
colMeans(rdt) |> sort()
# mean time
dt |> group_by(test) |> summarise(mean_time = mean(time))
