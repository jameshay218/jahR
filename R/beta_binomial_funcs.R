#' Beta binomial density
#'
#' @param x the number of successes
#' @inheritParams pbb
#' @return the probability density of the beta binomial with given input parameters
#' @family beta_binomial_functions
#' @examples
#' dbb(10,50,1,1)
#' @export
dbb <- function(x, u, v, N) {
  beta(x + u, N - x + v) / beta(u, v) * choose(N, x)
}

#' Beta binomial distribution function
#'
#' @param q vector of quantiles
#' @param N the number of trials
#' @param u equivalent to alpha in description
#' @param v equivalent to beta in description
#' @return the beta binomial distribution function
#' @family beta_binomial_functions
#' @examples
#' pbb(10,50,1,1)
#' @export
pbb <- function(q, u, v, N) {
  sapply(q, function(xx) sum(dbb(0:xx, N, u, v)))
}

#' Beta binomial quantile function
#'
#' @param p vector of probabilities
#' @inheritParams pbb
#' @return the beta binomial quantile function
#' @family beta_binomial_functions
#' @examples
#' qbb(10,50,1,1)
#' @export
qbb <- function(p, u, v,N) {
  pp <- cumsum(dbb(0:N, N, u, v))
  sapply(p, function(x) sum(pp < x))
}

#' Beta binomial random generation
#'
#' @param n number of draws
#' @inheritParams pbb
#' @return the beta binomial random draw function
#' @family beta_binomial_functions
#' @examples
#' rbb(100,50,1,1)
#' @export
rbb <- function(n, u, v, N) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}

#' Beta binomial mean
#'
#' @param n the number of trials
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @return calculates the mean of the beta binomial function
#' @family beta_binomial_functions
#' @examples
#' bb_mean(100, 1, 1)
#' @export
bb_mean <- function(n, alpha, beta) {
  return(n * alpha / (alpha + beta))
}

#' Beta binomial variance
#'
#' @inheritParams bb_mean
#' @return calculates the variance of the beta binomial function
#' @examples
#' bb_var(100, 1, 1)
#' @export
bb_var <- function(n, alpha, beta) {
  top <- n * alpha * beta * (alpha + beta + n)
  bot <- ((alpha + beta)^2) * (alpha + beta + 1)
  return(top / bot)
}



#' Beta binomial histogram
#'
#' Plots a histogram of the beta binomial function with given mean and variance
#' @inheritParams find_prior_alpha_beta
#' @family beta_binomial_functions
#' @examples
#' hist_rbb(100, 0.15, 0.1)
#' @export
hist_rbb <- function(n, mean, var) {
  pars <- find_prior_alpha_beta(n, mean, var)
  a <- pars["a"]
  b <- pars["b"]
  hist(rbb(10000, n, a, b), breaks = seq(-1, n, by = 1))
}

#' Beta binomial parameters
#'
#' Finds the required alpha and beta to give a desired mean and variance of the beta binomal function
#' @param n the number of trials
#' @param mean desired mean
#' @param var desired variance
#' @return alpha and beta for the beta binomial
#' @family beta_binomial_functions
#' @examples
#' find_prior_alpha_beta(n=100,mean=25,var=5)
#' @export
find_prior_alpha_beta <- function(n, mean, var) {
  y <- mean
  z <- var
  x <- (n - y) / y
  top <- z * (1 + x)^2 - (n^2) * x
  bot <- n * x * (1 + x) - z * ((1 + x)^3)
  a2 <- top / bot
  b2 <- x * a2
  return(c("a" = a2, "b" = b2))
}

#' Beta binomial parameter match
#'
#' Given a beta binomial function with known N, alpha and beta, finds the alpha and beta for a beta binomial with a different N that gives the same mean and variance
#' @param n1 the number of trials in the first distribution
#' @param a1 alpha in the first distribution
#' @param b1 beta in the first distribution
#' @param n2 the number of trials in the second distribution
#' @return alpha and beta for the second distribution
#' @family beta_binomial_functions
#' @export
find_bb_2 <- function(n1, a1, b1, n2) {
  y <- bb_mean(n1, a1, b1)
  z <- bb_var(n1, a1, b1)
  x <- (n2 - y) / y
  top <- z * (1 + x)^2 - (n2^2) * x
  bot <- n2 * x * (1 + x) - z * ((1 + x)^3)
  a2 <- top / bot
  b2 <- x * a2
  return(c("a2" = a2, "b2" = b2))
}


