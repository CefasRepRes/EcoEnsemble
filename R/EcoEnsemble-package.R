#' Bayesian ensemble functions for combining ecosystem models
#'
#' @description The `EcoEnsemble` package implements the general ensemble framework for
#' ecosystem models laid out in Spence et al (2018) in R.
#'
#' @details
#' The ensemble model can be implemented in three main stages:
#' 1. Eliciting priors on discrepancy terms: This is done by using the
#'  `define_priors` function.
#' 2. Fitting the ensemble model: Using `fit_ensemble_model` with simulator
#' outputs, observations, and priors, the ensemble model can be fit, obtaining
#' either the point estimate which maximises the posterior density, or running MCMC to
#' generate a sample from the posterior denisty of the ensemble model.
#' 3. Sampling from the Kalman filter: Using `generate_sample` with
#' the fitted ensemble object, the discrepancy terms and the ensemble's best
#' guess of the truth can be generated. Similarly to `fit_ensemble_model`, this
#' can either be a point estimate (the MLE) or a full sample of the Kalman filter.
#'
#' @docType package
#' @name EcoEnsemble-package
#' @aliases EcoEnsemble
#' @useDynLib EcoEnsemble, .registration = TRUE
#' @import methods
#' @import matrixcalc
#' @import Rcpp
#' @importFrom rstan sampling optimizing
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#' Spence et. al. (2018). A general framework for combining ecosystem models. https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310
NULL