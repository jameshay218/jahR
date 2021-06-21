#' Generalized logistic function
#'
#' Implementation of the fully generalized logistic function, as described at \url{https://en.wikipedia.org/wiki/Generalised_logistic_function}.
#' Note that the default arguments are set to give the logistic function, but each parameter can be customized.
#' @param x the vector of values to solve the function over
#' @param A the lower asymptote
#' @param K the upper asymptote, also the carrying capacity when A=0 and C=1
#' @param C Usually C=1, but if C != 1 then the upper asymptote = A + (K-A)/(C^(1/v))
#' @param Q related to the initial value at Y(0)
#' @param k IMPORTANT: the growth rate
#' @param M the starting time ie. shifts in the x-axis
#' @param v v>0, affecting where maximum growth occurs
#' @return a vector of values for each entry of x, giving y(x) where y() is the generalized logistic function.
#' @export
#' @examples
#' y <- generalized_logistic_function(x=seq(-10,10,by=0.1), k=0.5)
generalized_logistic_function <- function(x,A=0,K=1,C=1,Q=1,k=0,M=0,v=1){
  A + (K-A)/(C+Q*exp(-k*(x-M))^(1/v))
}
