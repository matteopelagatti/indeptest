#' Pelagatti-Monti test of idependence
#' 
#' It tests the idependences of two random variables x and y, which can be
#' numerical or categorical (factor).
#' 
#' @param x either a numerical vector or a factor.
#' @param y either a numerical vector or a factor (of the same length of x).
#' @param order integer or vector of two integers; if determines the number of basis
#' functions to approximate any f(x) and g(y); if integer both functions use the same
#' number of basis functions otherwise, the first integer determines the order for f(x)
#' and the second integer determines the order of g(y).
#' @param basis string or vector of two strings; the type of basis functions to be used
#' for f(x) and g(y); if only one string is passed the same basis is used both for f(x)
#' and g(y); the choices are ("spline", "orthopoly", "poly", "dummy"), but dummy is
#' available only for categorial variables (factor).
#' @param ties.method string among ("random", "average", "first", "max", "min") passed
#' to the ranking function; if random is used, the statistic may change value if ties
#' are present, however this assures that the p-value is correcly computed (uniformely
#' distributed). Do not use "first": very strong size distortion.
#' 
#' @return A list containing:\describe{
#'    \item{stat}{The test statistic}
#'    \item{pvalue}{The pvalue based on the asymptotic approximation}
#'    \item{n}{The size of the sample}
#'    \item{order}{A vector of integers with the number of basis functions
#'                 used for x and y}
#'    \item{basis}{A vector of strings with the names of the basis functions
#'                 used for x and y}
#'    \item{cor}{A vector with the square root of the eigenvalues: the first
#'               element is the sample Renyi maximal correlation.}
#'  }
#'  
#' @examples
#' indeptest(rnorm(100), rnorm(100), 3, "spline")
#' indeptest(rnorm(100), sample(c("yes", "no"), 100, TRUE), 3, c("spline", "dummy"))
indeptest <- function(x, y, order = max(1, floor((NROW(x))^(1/3)) - 1),
                      basis = c("spline", "orthopoly", "poly", "dummy"),
                      ties.method = c("random", "average", "first", "max", "min")) {
  ties.method <- match.arg(ties.method)
  basis <- match.arg(basis, several.ok = TRUE)
  n <- length(x)
  if (length(y) != n) stop("x and y must have equal length")
  if (is.character(x)) x <- factor(x)
  if (is.character(y)) y <- factor(y)
  if (length(basis) == 2) {
    basisx <- basis[1]
    basisy <- basis[2]
  } else {
    basisx <- basisy <- basis[1]
  }
  if (length(order) == 2) {
    orderx <- order[1]
    ordery <- order[2]
  } else {
    orderx <- ordery <- order[1]
  }
  if (is.factor(x) & basisx == "dummy" ) {
    levx <- unique(x)
    X <- outer(x, levx[-1], `==`)
    orderx <- length(levx) - 1
  } else {
    Rx <- rank(x, ties.method = ties.method)/(n + 1)
  }
  if (is.factor(y) & basisy == "dummy") {
    levy <- unique(y)
    Y <- outer(y, levy[-1], `==`)
    ordery <- length(levy) - 1
  } else {
    Ry <- rank(y, ties.method = ties.method)/(n + 1)
  }
  if (is.factor(x) && (basisx != "dummy") && (orderx >= length(unique(x)))) {
    orderx <- length(unique(x)) - 1
    warning("The variable x is categorical and the number of basis functions must 
         be smaller than the number of categories: order_x was set to ", orderx)
  }
  if (is.factor(y) && (basisy != "dummy") && (ordery >= length(unique(y)))) {
    ordery <- length(unique(y)) - 1
    warning("The variable y is categorical and the number of basis functions must 
         be smaller than the number of categories: order_y was set to ", ordery)
  }
  if (basisx == "spline") {
    degreex <- if (orderx < 3) orderx else 3
    X <- splines::bs(Rx, df = orderx, degree = degreex)
  }
  if (basisy == "spline") {
    degreey <- if (ordery < 3) ordery else 3
    Y <- splines::bs(Ry, df = ordery, degree = degreey)
  }
  if (basisx == "orthopoly") X <- poly(Rx, orderx)
  if (basisy == "orthopoly") Y <- poly(Ry, ordery)
  if (basisx == "poly")  X <- outer(Rx, 1:orderx, "^")
  if (basisy == "poly")  Y <- outer(Ry, 1:ordery, "^")
  # --- using cancor(), slowest ---
  # cc <- cancor(X, Y)
  # if (order != length(cc$cor)) {
  #   warning("Unable to compute basis of order ", order, ". Order was set to ", length(cc$cor))
  #   order <- length(cc$cor)
  # }
  # stat <- -(n - order - 1.5)*sum(log(1 - cc$cor^2))
  # --- As cancor() but without overhead and eigenvectors, slow ---
  # XX <- X - rep(colMeans(X), rep.int(n, order))
  # YY <- Y - rep(colMeans(Y), rep.int(n, order))
  # qx <- qr(XX)
  # qy <- qr(YY)
  # dx <- qx$rank
  # dy <- qy$rank
  # cc <- svd(qr.qty(qx, qr.qy(qy, diag(1, n, dy)))[1L:dx, , drop = FALSE], dx, dy)$d
  # if (order != length(cc)) {
  #   warning("Unable to compute basis of order ", order, ". Order was set to ", length(cc))
  #   order <- length(cc)
  # }
  # stat <- -(n - order - 1.5)*sum(log(1 - cc^2))
  # --- using eigenvalues of product of cov matrices, faster ---
  # Sxx <- var(X)
  # Sxy <- cov(X, Y)
  # Syy <- var(Y)
  # cc2 <- eigen(chol2inv(chol(Sxx, pivot = T)) %*% Sxy %*% chol2inv(chol(Syy, pivot = T)) %*% t(Sxy),
  #         symmetric = FALSE, only.values = TRUE)$values
  # cc <- sqrt(cc2)
  # stat <- -(n - order - 1.5)*sum(log(1 - cc2))
  # --- similar to above, but should be more stable, fast ---
  # Sxx <- var(X)
  # Sxy <- cov(X, Y)
  # Syy <- var(Y)
  # e <- eigen(Sxx)
  # invSqrtSxx <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)
  # H <- invSqrtSxx %*% Sxy
  # cc2 <- eigen(H %*% chol2inv(chol(Syy)) %*% t(H), symmetric = TRUE, only.values = TRUE)$values
  # cc <- sqrt(cc2)
  # stat <- -(n - order - 1.5)*sum(log(1 - cc2))
  # --- like above, but using a generalized eigenvalues algorithm
  Sxx <- var(X)
  Sxy <- cov(X, Y)
  Syy <- var(Y)
  min_order <- min(orderx, ordery)
  max_order <- max(orderx, ordery)
  if (dim(Sxx)[1] == 1) {
    cc2 <- Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy) / Sxx
  } else {
    cc2 <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx, symmetric = T, 
                          only.values = TRUE)$values[(max_order - min_order + 1):(max_order)]
  }
  cc2[cc2<0] <- 0
  cc <- sqrt(sort(cc2, T))
  stat <- - (n - (orderx + ordery + 3)/2) * sum(log(1 - cc2))
  
  list(stat = stat, pvalue = pchisq(stat, orderx * ordery, lower.tail = FALSE), 
       n = n, order = c(order_x = orderx, order_y = ordery),
       basis = c(basis_x = basisx, basis_y = basisy), cor = cc)
}

cat_test <- function(n, ncatx, ncaty, order = max(1, floor((NROW(x))^(1/3)) - 1),
                     basis = c("spline", "orthopoly", "poly", "dummy"),
                     ties.method = c("random", "average", "first", "max", "min")) {
  x <- sample(as.character(1:ncatx), n, TRUE)
  y <- sample(as.character(1:ncaty), n, TRUE)
  tst <- indeptest(x, y, order, basis, ties.method)
  c(stat = tst$stat, pvalue = tst$pvalue, tst$order[1], tst$order[2])
}

sim_cat_dist <- function(nsim, n, ncatx, ncaty, order, basis, ties.method) {
  sims <- replicate(nsim, cat_test(n, ncatx, ncaty, order, basis, ties.method))
  hist(sims["stat", ], freq = FALSE)
  curve(dchisq(x, df = sims["order_x", 1]*sims["order_y", 1]), add = TRUE)
  hist(sims["pvalue", ], freq = FALSE)
  lines(x = c(0, 1),  y = c(1, 1))
  invisible(sims)
}

depanalysis <- function(x, y, order = max(1, floor((NROW(x))^(1/3)) - 1),
                        basis = c("spline", "orthopoly", "poly", "dummy"),
                        ties.method = c("random", "average", "first", "max", "min"),
                        dimensions = 1:(min(4, order))) {
  ties.method <- match.arg(ties.method)
  basis <- match.arg(basis, several.ok = TRUE)
  n <- length(x)
  if (length(y) != n) stop("x and y must have equal length")
  if (is.character(x)) x <- factor(x)
  if (is.character(y)) y <- factor(y)
  if (length(basis) == 2) {
    basisx <- basis[1]
    basisy <- basis[2]
  } else {
    basisx <- basisy <- basis[1]
  }
  if (length(order) == 2) {
    orderx <- order[1]
    ordery <- order[2]
  } else {
    orderx <- ordery <- order[1]
  }
  if (is.factor(x) & basisx == "dummy" ) {
    levx <- unique(x)
    X <- outer(x, levx[-1], `==`)
    orderx <- length(levx) - 1
  } else {
    Rx <- rank(x, ties.method = ties.method)/(n + 1)
  }
  if (is.factor(y) & basisy == "dummy") {
    levy <- unique(y)
    Y <- outer(y, levy[-1], `==`)
    ordery <- length(levy) - 1
  } else {
    Ry <- rank(y, ties.method = ties.method)/(n + 1)
  }
  if (is.factor(x) && (basisx != "dummy") && (orderx >= length(unique(x)))) {
    orderx <- length(unique(x)) - 1
    warning("The variable x is categorical and the number of basis functions must 
         be smaller than the number of categories: order_x was set to ", orderx)
  }
  if (is.factor(y) && (basisy != "dummy") && (ordery >= length(unique(y)))) {
    ordery <- length(unique(y)) - 1
    warning("The variable y is categorical and the number of basis functions must 
         be smaller than the number of categories: order_y was set to ", ordery)
  }
  if (basisx == "spline") {
    degreex <- if (orderx < 3) orderx else 3
    X <- splines::bs(Rx, df = orderx, degree = degreex)
  }
  if (basisy == "spline") {
    degreey <- if (ordery < 3) ordery else 3
    Y <- splines::bs(Ry, df = ordery, degree = degreey)
  }
  if (basisx == "orthopoly") X <- poly(Rx, orderx)
  if (basisy == "orthopoly") Y <- poly(Ry, ordery)
  if (basisx == "poly")  X <- outer(Rx, 1:orderx, "^")
  if (basisy == "poly")  Y <- outer(Ry, 1:ordery, "^")

  canc <- cancor(X, Y)
  cc2 <- canc$cor^2
  cc2[cc2<0] <- 0
  stat <- - (n - (orderx + ordery + 3)/2) * sum(log(1 - cc2))

  # plots
  ordx <- order(x)
  ordy <- order(y)
  fx <- as.matrix(X %*% canc$xcoef[, dimensions])
  fy <- as.matrix(Y %*% canc$ycoef[, dimensions])
  colnames(fx) <- paste0("dim", dimensions)
  colnames(fy) <- paste0("dim", dimensions)
  
  dfx <- data.frame(x = x[ordx], fx[ordx, ])
  dfy <- data.frame(y = y[ordy], fy[ordy, ])
  dfxy <- cbind(tidyr::pivot_longer(data.frame(fx), everything(),
                                    names_to = "dim", values_to = "f(x)"),
                tidyr::pivot_longer(data.frame(fy), everything(),
                                    names_to = "dimy", values_to = "f(y)")
  )
  dfxy$dim <- as.factor(dfxy$dim)
  levels(dfxy$dim) <- paste(levels(dfxy$dim),
                            " cor =",
                            round(canc$cor[dimensions], 2))
  names(dfx) <- c("x", paste0("dim", dimensions))
  names(dfy) <- c("y", paste0("dim", dimensions))

  dfx <- tidyr::pivot_longer(dfx, starts_with("dim"),
                             names_to = "dim",
                             values_to = "f(x)")
  dfy <- tidyr::pivot_longer(dfy, starts_with("dim"),
                             names_to = "dim",
                             values_to = "f(y)")
  
  px <- ggplot2::ggplot(dfx, ggplot2::aes(x = x, y = `f(x)`)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(dim), scales = "free_y")
  py <- ggplot2::ggplot(dfy, ggplot2::aes(x = y, y = `f(y)`)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(dim), scales = "free_y")

  pxy <- ggplot2::ggplot(dfxy, ggplot2::aes(x = `f(x)`, y = `f(y)`)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = "y~x", se = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(dim), scales = "free")

  list(stat = stat, pvalue = pchisq(stat, orderx * ordery, lower.tail = FALSE), 
       n = n, order = c(order_x = orderx, order_y = ordery),
       basis = c(basis_x = basisx, basis_y = basisy), cor = canc$cor,
       plotx = px, ploty = py, plotxy = pxy)
}