### Examples from 
### Liu et al. (2022) Model-free Feature Screening and FDR Control with Knockoff Features



# Example 1 ---------------------------------------------------------------
## Linear models
## for each model the active set contains the first 5 covariates in x.

# Model 1.a: x \sim N(0;Sigma) and varepsilon \sim N(0; 1).
# Model 1.b: x \sim N(0;Sigma) and varepsilon \sim Cauchy(0; 1).
### n.b. Cauchy \sim (0,1) \sim t(df =1) Student's t distribution

# Model 1.c: u \sim Cauchy(0; Ip), x = Sigma^{1/2}u and varepsilon \sim N(0; 1).
# Model 1.d: u \sim Cauchy(0; Ip), x = Sigma^{1/2}u and varepsilon \sim Cauchy(0; 1).
# Model 1.e: (Continuous) Y = exp(x^T \beta)+ varepsilon, where varepsilon \sim N(0; 1).
# Model 1.f: (Discrete) Y \sim Poisson(exp(x^T \beta)).

require(MASS)
require(mvtnorm)

## CASO 1 
# Model 1.a e 1.b in Liu et al.(2022).
GendataLM <- function (n, p, rho, beta = c(rep(1, 5), rep(0, p - 5)), 
                       error = c("gaussian", "t", "cauchy")) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  if (error == "gaussian" | is.null(error)) {
    Y = X %*% beta + rnorm(n)
  }
  else if (error == "t") {
    Y = X %*% beta + rt(n, 2)
  }
  else if (error == "cauchy") {
    Y = X %*% beta + rcauchy(n)
  }
  else {
    stop("The author has not implemented this error term yet.")
  }
  return(list(X = X, Y = Y))
}

## CASO 2
# Model 1.c 1.d in Liu et al.(2022).
# heavy-tailed
GendataHT <- function (n, p, rho, beta = c(rep(1, 5), rep(0, p - 5)),
                       error = c("gaussian", "cauchy")) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  U <- rmvt(n, sigma = diag(p), df = 1)
  X = U %*% sqrt(sig)
  if (error == "gaussian") {
    Y = X %*% beta + rnorm(n)
  } else if (error == "cauchy") {
    Y = X %*% beta + rcauchy(n)
  }
  else {
    stop("The author has not implemented this error term yet.")
  }
  return(list(X = X, Y = Y))
}
## GendataHT(10,50, 0.5, error = "gaussian")

# Model 1.e in Liu et al.(2022).
## CASO  3

Gendata1e <- function (n, p, rho, beta = c(rep(2, 5), rep(0, p - 5))) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = exp(X %*% beta) + rnorm(n)
  return(list(X = X, Y = Y))
}


# Model 1.f in Liu et al.(2022).
# Usage GendataPM(n, p, rho, beta = c(rep(1, 5), rep(0, p - 5)))

## CASO  4

GendataPM <- function (n, p, rho, beta = c(rep(2, 5), rep(0, p - 5))) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  myrates = exp(X %*% beta)
  Y = rpois(n, myrates)
  return(list(X = X, Y = Y))
}


# Example 2 ---------------------------------------------------------------
## NonLinear models
## for each model the active set contains the first 4 covariates in x.

# Model 2.a: Y = 5*X1 + 2*sin(pi*X2/2) + 2*X3*I(X3 > 0) + 2*exp(5*X4)
#                + varepsilon with varepsilon \sim N(0; 1).

# Model 2.b: Y = 3*X1 + 3*X2^3 + 3X3^(-1) + 5*I(X4 > 0) 
#                + varepsilon with varepsilon \sim N(0; 1).

# Model 2.c: Y = 1 - 5(X2 + X3)^3 * exp(-5(X1 + X4^2 )) 
#                + varepsilon with varepsilon \sim N(0; 1).

# Model 2.d: Y = 1 - 5(X2 + X3)^(-3)* exp(1 + 10 sin( \pi X1/2) + 5X4) 
#                + varepsilon with varepsilon \sim N(0; 1).

## Model 2.e: Y = x^T \beta exp(x^T \beta)+ varepsilon, where varepsilon \sim N(0; 1).


Gendata2a <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 5*X[ ,1] + 2*sin(pi*X[ ,2]/2) + 2*X[ ,3]*(X[ ,3] > 0) + 2*exp(5*X[ ,4]) + rnorm(n)
  return(list(X = X, Y = Y))
}

## Gendata2a(n=10, p=50, rho=0.5)

Gendata2b <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 3*X[ ,1] + 3*X[ ,2]^3 + 3*1/X[ ,3] + 5*(X[ ,4] > 0)  + rnorm(n)
  return(list(X = X, Y = Y))
}

##  Gendata2b(n=10, p=50, rho=0.5)

Gendata2c <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 1 - 5*(X[ ,2] + X[ ,3])^3 * exp(-5*(X[ ,1] + X[ ,4]^2 ))   + rnorm(n)
  return(list(X = X, Y = Y))
}

##  Gendata2c(n=10, p=50, rho=0.5)

Gendata2d <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 1 - 5*(X[ ,2] + X[ ,3])^(-3) *exp(1 + 10*sin(pi*X[ ,1]/2) + 5*X[ ,4]) + rnorm(n)
  return(list(X = X, Y = Y))
}

##  Gendata2d(n=10, p=50, rho=0.5)

## multiple index model (Zhu et al. 2011) Example 3.b 

## CASO 6 

Gendata2e <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  u1 <- runif(1)
  u2 <- runif(1)
  u3 <- runif(1)
  u4 <- runif(1)
  beta1 = c(2-u1,2-u2, rep(0, p - 2))
  beta2 = c(0,0, 2+u3, 2+u4, rep(0, p-4))
  Y = X %*% beta1 + exp(X %*% beta2) + rnorm(n)
  return(list(X = X, Y = Y))
}


## heteroschedastic model (Zhu et al. 2011) (Zhu et al. 2011) Example 3.c

## CASO 7

Gendata2f <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  u1 <- runif(1)
  u2 <- runif(1)
  u3 <- runif(1)
  u4 <- runif(1)
  beta1 = c(2-u1,2-u2, rep(0, p - 2))
  beta2 = c(0,0, 2 + u3, 2 + u4, rep(0, p - 4))
  Y = X %*% beta1 + exp(X %*% beta2 + rnorm(n)) 
  return(list(X = X, Y = Y))
}


## transformation model (Zhu et al. 2011) (Zhu et al. 2011) similar to Example 3.a

## CASO 5

Gendata2g <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  u1 <- runif(1)
  u2 <- runif(1)
  u3 <- runif(1)
  u4 <- runif(1)
  beta = c(2 - u1, 2 - u2, 2 - u3, 2 - u4 ,rep(0, p - 4))
  Y = exp(X %*% beta + rnorm(n)) 
  return(list(X = X, Y = Y))
}


##  Case 8

Gendata2c <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 1 - 5*(X[ ,2] + X[ ,3])^3 * exp(-5*(X[ ,1] + X[ ,4]^2 ))   + rnorm(n)
  return(list(X = X, Y = Y))
}

##  Case 9

Gendata2d <- function (n, p, rho) {
  sig = matrix(0, p, p)
  sig = rho^abs(row(sig) - col(sig))
  diag(sig) <- rep(1, p)
  X = mvrnorm(n, rep(0, p), sig)
  Y = 1 - 5*(X[ ,2] + X[ ,3])^(-3) *exp(1 + 10*sin(pi*X[ ,1]/2) + 5*X[ ,4]) + rnorm(n)
  return(list(X = X, Y = Y))
}

