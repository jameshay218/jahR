% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antibody_kinetics.R
\name{plot_model_fits_mean}
\alias{plot_model_fits_mean}
\title{Plots posterior estimates for the population-mean antibody trajectory}
\usage{
plot_model_fits_mean(
  model_func,
  chain,
  parTab,
  nsamp = 100,
  times = seq(0, 365, by = 1),
  use_log = TRUE,
  plot_draws = TRUE
)
}
\arguments{
\item{model_func}{pointer to antibody kinetics function used to general modeled trajectories}

\item{chain}{the MCMC chain object to draw posterior samples from}

\item{times}{vector of times to solve model over if expand_times is TRUE}

\item{use_log}{if TRUE, log10 transforms the model-predicted titers}

\item{plot_draws}{if TRUE, then plots individual draws from the posterior, otherwise plots a ribbon of credible/prediction intervals and posterior median.}
}
\value{
a ggplot2 object
}
\description{
Plots posterior estimates for the population-mean antibody trajectory
}
