#' A general framework for combining ecosystem models
#'
#' @description The `EcoEnsemble` package implements the framework for combining ecosystem models laid out in Spence et al (2018).
#'
#' @details
#' The ensemble model can be implemented in three main stages:
#' 1. Eliciting priors on discrepancy terms: This is done by using the `EnsemblePrior` constructor.
#' 2. Fitting the ensemble model: Using `fit_ensemble_model` with simulator outputs, observations and prior information. The ensemble model can be fit, obtaining either the point estimate, which maximises the posterior density, or running Markov chain Monte Carlo to generate a sample from the posterior denisty of the ensemble model.
#' 3. Sampling the latent variables from the fitted model: Using `generate_sample` with the fitted ensemble object, the discrepancy terms and the ensemble's best guess of the truth can be generated. Similarly to `fit_ensemble_model`, this can either be a point estimate or a full sample.
#'
#' @useDynLib EcoEnsemble, .registration = TRUE
#' @import Rcpp
#' @import methods
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import tibble
#' @import cowplot
#' @import matrixcalc
#' @importFrom rstan sampling optimizing
#'
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. [https://mc-stan.org](https://mc-stan.org)
#' @references Spence, M. A., J. L. Blanchard, A. G. Rossberg, M. R. Heath, J. J. Heymans, S. Mackinson, N. Serpetti, D. C. Speirs, R. B. Thorpe, and P. G. Blackwell. 2018. "A General Framework for Combining Ecosystem Models." Fish and Fisheries 19: 1013-42. [https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310](https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310)
"_PACKAGE"
