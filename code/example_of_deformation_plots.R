library(splines)
library(ggplot2)
source("indeptest_with_cat.R")

n <- 1000
set.seed(2022)

# Example 1: x is the standard deviation of y
x1 <- rchisq(n, df = 2) |> sqrt()
y1 <- rnorm(n)*x1

ex1_raw_plot <- ggplot(data.frame(x = x1, y = y1,
                                  lab = paste("cor = ", round(cor(x1, y1), 2))),
                       aes(x = x, y = y)) +
                  geom_point() +
                  geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
                  facet_wrap(~lab) +
                  ggtitle("Example 1: scatter of raw data")

m1 <- depanalysis(x1, y1, basis = "orthopoly", dimensions = 1)

m1$plotx + ggtitle("Example 1: x-transforms")
m1$ploty + ggtitle("Example 1: y-transforms")
m1$plotxy + ggtitle("Example 1: scatter of the transforms")

# Example 2: circle with noise
theta <- runif(n, 0, 2*pi)
x2 <- cos(theta) + rnorm(n, sd = 0.1)
y2 <- sin(theta) + rnorm(n, sd = 0.1)

ex2_raw_plot <- ggplot(data.frame(x = x2, y = y2,
                                  lab = paste("cor = ", round(cor(x2, y2), 2))),
                       aes(x = x, y = y)) +
                  geom_point() +
                  geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
                  facet_wrap(~lab) +
                  ggtitle("Example 2: scatter of raw data")

ggplot(data.frame(x = x2, y = y2), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
  ggtitle("Example 2: raw data")

m2 <- depanalysis(x2, y2, basis = "orthopoly", dimensions = 1)

m2$plotx + ggtitle("Example 2: x-transforms")
m2$ploty + ggtitle("Example 2: y-transforms")
m2$plotxy + ggtitle("Example 2: scatter of the transforms")

# Example 3: linear dependence only on a small intervall
x3 <- rnorm(n)
y3 <- ifelse(x3 < -2, -x3, ifelse(x3 > 2, x3, runif(n, -2, 2)))

ex3_raw_plot <- ggplot(data.frame(x = x3, y = y3,
                                  lab = paste("cor = ", round(cor(x3, y3), 2))),
                       aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
  facet_wrap(~lab) +
  ggtitle("Example 3: scatter of raw data")

m3 <- depanalysis(x3, y3, order = 15, basis = "spline", dimensions = 1)

m3$plotx + ggtitle("Example 3: x-transforms")
m3$ploty + ggtitle("Example 3: y-transforms")
m3$plotxy + ggtitle("Example 3: scatter of the transforms")

w <- 4
h <- 3
# example 1
ggsave("ex1_raw.pdf",
       ex1_raw_plot +
         theme_bw() +
         ggtitle("Example 1: scatter of raw data"),
       width = w, height = h)
ggsave("ex1_xtrans.pdf",
       m1$plotx +
         theme_bw() +
         ylab("g(x)") +
         ggtitle("Example 1: x-transform"),
       width = w, height = h)
ggsave("ex1_ytrans.pdf",
       m1$ploty +
         theme_bw() +
         ylab("h(y)") +
         ggtitle("Example 1: y-transform"),
       width = w, height = h)
ggsave("ex1_fxfy.pdf",
       m1$plotxy +
         theme_bw() +
         xlab("g(x)") + ylab("h(y)") +
         ggtitle("Example 1: scatter of the transforms"),
       width = w, height = h)

# example 2
ggsave("ex2_raw.pdf",
       ex2_raw_plot +
         theme_bw() +
         ggtitle("Example 2: scatter of raw data"),
       width = w, height = h)
ggsave("ex2_xtrans.pdf",
       m2$plotx +
         theme_bw() +
         ylab("g(x)") +
         ggtitle("Example 2: x-transform"),
       width = w, height = h)
ggsave("ex2_ytrans.pdf",
       m2$ploty +
         theme_bw() +
         ylab("h(y)") +
         ggtitle("Example 2: y-transform"),
       width = w, height = h)
ggsave("ex2_fxfy.pdf",
       m2$plotxy +
         theme_bw() +
         xlab("g(x)") + ylab("h(y)") +
         ggtitle("Example 2: scatter of the transforms"),
       width = w, height = h)

# example 3
ggsave("ex3_raw.pdf",
       ex3_raw_plot +
         theme_bw() +
         ggtitle("Example 3: scatter of raw data"),
       width = w, height = h)
ggsave("ex3_xtrans.pdf",
       m3$plotx +
         theme_bw() +
         ylab("g(x)") +
         ggtitle("Example 3: x-transform"),
       width = w, height = h)
ggsave("ex3_ytrans.pdf",
       m3$ploty +
         theme_bw() +
         ylab("h(y)") +
         ggtitle("Example 3: y-transform"),
       width = w, height = h)
ggsave("ex3_fxfy.pdf",
       m3$plotxy +
         theme_bw() +
         xlab("g(x)") + ylab("h(y)") +
         ggtitle("Example 3: scatter of the transforms"),
       width = w, height = h)
