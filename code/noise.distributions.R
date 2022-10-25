# the X values were simulated iid from a uniform distribution
# Linear: f (x) = x , error distribution is N(0; 1)
# Quadratic: f (x) = 4(x - 1/2)^2 , error distribution is N(0; 1)
# Cubic: f (x) = 128(x - 1/3)^3 - 48(x - 1/3)^2 - 12(x - 1/3), error distribution is N(0; 100)
# Sine: f (x) = sin(4*pi*x) , error distribution is N(0; 4)
# X^1/4: f (x) = x*1/4 , error distribution is N(0; 1)
# Circle: f (x) = (2r - 1) sqrt(1 -(2x- 1)^2) , error distribution is N(0; 1/16), where r is a Bernoulli(1/2) variable
# Two curves: f (x) = 2*r*x + (1 - r )*sqrt(x)/2 , error distribution is N(0; 1/4), where r is a Bernoulli(1/2) variable
# X-function: f (x) = r*x+(1-r)*(1-x), error distribution is N(0; 1/25), where r is a Bernoulli(1/2) variable
# Diamond: f (x) = r1 * (x < 0.5) + r2 * (x >= 0.5), error distribution is N(0; 1/100), where r1 is a U(0.5- x; 0.5 + x)
# variable and r2 is a U(x - 0.5; 1.5- x) variable

datagen.noise = function(n, typ, noise=.1) {
  
  # lin+noise          typ==1
  # parabolic+noise    typ==2
  # cubic+noise        typ==3
  # sin+noise          typ==4
  # x^(1/4) + noise    typ==5
  # circle             typ==6
  # two curves         typ==7
  # X function         typ==8
  # Diamond            typ==9
  # XsdY               typ==10
  
  if (typ == '') {
  } else if (typ == 1) {
    linear(n, noise)
  } else if (typ == 2) {
    quadratic(n, noise)
  } else if (typ == 3) {
    cubic(n, noise)
  } else if (typ == 4) {
    sine(n, noise)
  } else if (typ == 5) {
    x14(n, noise)
  } else if (typ == 6) {
    circle(n, noise)
  } else if (typ == 7) {
    twocurves(n, noise)
  } else if (typ == 8) {
    Xfun(n, noise)
  } else if (typ == 9) {
    Diamond(n, noise)
  } else if (typ == 10) {
    XsdY(n, noise)
  }
}

linear <- function(n, noise){ # 1
  x <- runif(n)
  y <- x + noise * rnorm(n)
  return(rbind(x, y))
}

quadratic <- function(n, noise){ # 2
  x <- runif(n)
  y <- 4*(x - 1/2)^2 + noise * rnorm(n)
  return(rbind(x, y))
}

cubic <- function(n, noise){ # 3
  x <- runif(n)
  y <- 128*(x - 1/3)^3 - 48*(x - 1/3)^2 - 12*(x - 1/3) + noise * rnorm(n, sd = 3)
  return(rbind(x, y))
}

sine <- function(n, noise){ # 4
  x <- runif(n)
  y <- sin(4*pi*x) + noise * rnorm(n, sd =2)
  return(rbind(x, y))
}

x14 <- function(n, noise){ # 5
  x <- runif(n)
  y <- x^(1/4) + noise * rnorm(n, sd = 0.5)
  return(rbind(x, y))
}

circle <- function(n, noise){ # 6
  x <- runif(n)
  r <- rbinom(n, size = 1, prob = 0.5)
  y <- (2*r - 1) * sqrt(1 -(2*x- 1)^2) + noise * rnorm(n, sd =1/2)
  return(rbind(x, y))
}

twocurves <- function(n, noise){ # 7
  x <- runif(n)
  r <- rbinom(n, size = 1, prob = 0.5)
  y <- 2*r*x + (1 - r )*sqrt(x)/2 + noise * rnorm(n)
  return(rbind(x, y))
}

Xfun <- function(n, noise){#8
  x <- runif(n)
  r=rbinom(n,1,0.5)
  y=r*x+(1-r)*(1-x) + noise*rnorm(n, sd = 1/5)
  return(rbind(x, y))
}

Diamond <- function(n, noise){ # 9
  x  <- runif(n)
  r1 <- runif(n, 0.5-x, 0.5+x)
  r2 <- runif(n, x-0.5, 1.5-x)
  y <- r1 * (x < 0.5) + r2 * (x >= 0.5) + noise *rnorm(n, sd = 1/10)
  return(rbind(x, y))
}

XsdY <- function(n, noise){ # 10
  x <- sqrt(rchisq(n, 1))
  y <- rnorm(n)*x + rnorm(n)*noise
  return(rbind(x, y))
}
