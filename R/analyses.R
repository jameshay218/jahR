#' Median and quantiles
#'
#' Calculates the median and 95% quantiles from a vector
#' @param x the vector
#' @return vector with 95% quantiles and median
#' @export
median_quantile <- function(x){
    out <- quantile(x, probs = c(0.025,0.5,0.975))
    names(out) <- c("ymin","y","ymax")
    return(out)
}


#' Find beta distribution mode
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @return the mode of the distribution
#' @export
get_beta_mode <- function(alpha, beta){
  mode_par <- (alpha-1)/(alpha+beta-2)
  mode_par
}
