% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antibody_kinetics.R
\name{extract_individual_parameters}
\alias{extract_individual_parameters}
\title{Extracts posterior draws for the individual-level parameters}
\usage{
extract_individual_parameters(chain, parTab, nsamp)
}
\arguments{
\item{chain}{the MCMC chain object to draw posterior samples from}
}
\value{
a tibble with posterior draws
}
\description{
Extracts posterior draws for the individual-level parameters
}
