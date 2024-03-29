% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antibody_kinetics.R
\name{kinetics_gamma}
\alias{kinetics_gamma}
\title{Gamma function kinetics}
\usage{
kinetics_gamma(times, pars)
}
\arguments{
\item{times}{vector of times to solve model over}

\item{pars}{vector of named model parameters}
}
\value{
vector of antibody levels
}
\description{
Very simple scaled gamma function to represent kinetics. Used in Zhao et al. 2018, but not particular interpretable otherwise.
}
\details{
Model parameters:
\itemize{
\item{"y0"}{initial antibody level}
\item{"peak"}{parameter to scale function in y-axis}
\item{"shape"}{gamma distribution shape}
\item{"scale"}{gamma distribution scale}
}
}
\examples{
pars <- c("y0"=0,"peak"=300,"shape"=1.2,"scale"=50)
times <- seq(0,1000,by=0.1)
y <- kinetics_gamma(times,pars)
plot(times, y, type='l')
}
\seealso{
Other kinetics_models: 
\code{\link{kinetics_ode}()},
\code{\link{kinetics_piecewise}()},
\code{\link{kinetics_power_function}()},
\code{\link{kinetics_ranjeva}()},
\code{\link{kinetics_simple}()}
}
\concept{kinetics_models}
