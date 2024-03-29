# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Fast implementation of the power-law waning antibody kinetics function
#' @export
#' @inheritParams kinetics_power_function
kinetics_power_function_cpp <- function(times, pars) {
    .Call('_jahR_kinetics_power_function_cpp', PACKAGE = 'jahR', times, pars)
}

#' Fast observation error function continuous
#' Calculate the probability of a set of observed titres given a corresponding set of predicted titres assuming continuous, bounded observations. FAST IMPLEMENTATION
#' @name Fast observation error function continuous
#' @param theta NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable titre
#' @param obs NumericVector, the vector of observed log titres
#' @param predicted_titres NumericVector, the vector of predicted log titres
#' @param a vector of same length as the input data giving the probability of observing each observation given the predictions
#' @return a likelihood for each observed titre
#' @export
likelihood_func_fast_continuous <- function(theta, obs, predicted_titres) {
    .Call('_jahR_likelihood_func_fast_continuous', PACKAGE = 'jahR', theta, obs, predicted_titres)
}

