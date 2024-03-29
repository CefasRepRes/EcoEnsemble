% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnsembleSample-class.R
\name{EnsembleSample}
\alias{EnsembleSample}
\title{A constructor for the \code{EnsembleSample} object}
\usage{
EnsembleSample(ensemble_fit, mle, samples)
}
\arguments{
\item{ensemble_fit}{An \code{EnsembleFit} object containing the fitted ensemble model.}

\item{mle}{An \code{array} of dimension \eqn{T \times (M + 2)\times N_{sample}} containing MLE point estimates from the \code{ensemble_fit} object, where \eqn{T} is the total time, \eqn{M} is the number of simulators and \eqn{N_{sample}} is the number of samples. For each time step, the \code{t}th element of the array is a \code{matrix} where each column is a sample and the rows are the variables:
\deqn{\left( y^{(t)}, \eta^{(t)}, z_1^{(t)}, z_2^{(t)}, \ldots, z_M^{(t)}\right)'}
where \eqn{y^{(t)}} is the ensemble model's prediction of the latent truth value at time \eqn{t},
\eqn{\eta^{(t)}} is the shared short-term discrepancy at time \eqn{t},
\eqn{z_i^{(t)}} is the individual short-term discrepancy of simulator \eqn{i} at time \eqn{t}.}

\item{samples}{An \code{array} of dimension \eqn{T \times (M + 2)\times N_{sample}} containing samples from the \code{ensemble_fit} object, where \eqn{T} is the total time, \eqn{M} is the number of simulators and \eqn{N_{sample}} is the number of samples. For each time step, the \code{t}th element of the array is a \code{matrix} where each column is a sample and the rows are the variables:
\deqn{\left( y^{(t)}, \eta^{(t)}, z_1^{(t)}, z_2^{(t)}, \ldots, z_M^{(t)}\right)'}
where \eqn{y^{(t)}} is the ensemble model's prediction of the latent truth value at time \eqn{t},
\eqn{\eta^{(t)}} is the shared short-term discrepancy at time \eqn{t},
\eqn{z_i^{(t)}} is the individual short-term discrepancy of simulator \eqn{i} at time \eqn{t}.}
}
\value{
An object of class \code{EnsembleSample}
}
\description{
A constructor for the \code{EnsembleSample} class. These objects are generated automatically using the \code{generate_sample} function.
}
\seealso{
\code{\link{EnsembleSample}}, \code{\link{generate_sample}}
}
