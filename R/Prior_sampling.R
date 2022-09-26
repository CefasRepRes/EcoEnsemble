#'Generate samples of parameters from prior distribution
#'
#'Methods to generates samples of the parameters from the prior distribution of the ensemble model.
#'
#'@param priors An `EnsemblePrior` object specifying the prior distributions for the ensemble.
#'@param M A `numeric` that represents the number of simulators. The default is 1.
#'@param full_sample A `logical` that runs a full sampling of the prior density of the ensemble model if `TRUE`. If `FALSE`, returns the point estimate which maximises the prior density of the ensemble model.
#'@param control If creating a full sample, this is a named `list` of paramaters to control Stan's sampling behaviour. See the documentation of the `stan()` function in the `rstan` package for details. The default value is `list(adapt_delta = 0.95)`. If optimizing, this value is ignored.
#'@param ... Additional arguments passed to the function \code{rstan::sampling} or  \code{rstan::optimizing}.
#'@return A `list` containing two items named `samples` and `point_estimate`. If `full_sample==TRUE`, `samples` is a `stanfit` and `point_estimate` is a `NULL` object, else `samples` is a `NULL` and `point_estimate` is a `list` object.
#'@references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'@seealso \code{\linkS4class{EnsembleFit}}
#'@export
#'@examples
#'\donttest{
#' priors <- EnsemblePrior(4)
#' prior_density <- prior_ensemble_model(priors, M = 4)
#' }
prior_ensemble_model <- function(priors, M = 1,
                                    full_sample = TRUE, control = list(adapt_delta = 0.95), ...){
  stan_input <- priors@priors_stan_input
  stan_input$M <- M
  stan_input$N <- priors@d

  #Using hierarchical priors uses a different model. This speeds up the sampling enormously
  mod <- stanmodels$ensemble_prior
  if(stan_input$form_prior_ind_st == 3){
    mod <- stanmodels$ensemble_prior_hierarchical
  }

  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(mod, data=stan_input, control = control, ...)
  }else{
    point_estimate <- rstan::optimizing(mod, data=stan_input,as_vector=FALSE, ...)
  }
  return(list(samples=samples, point_estimate=point_estimate))
}




#'Generate samples of latent variables from prior predictive distribution
#'
#'Methods to generates samples of the latent variables from the prior predictive distribution of the ensemble model.
#'
#'@inheritParams EnsembleData
#'@param full_sample A `logical` that runs a full sampling of the prior density of the ensemble model if `TRUE`. If `FALSE`, returns the point estimate which maximises the prior density of the ensemble model.
#'@param sam_priors A `list` containing two items named `samples` and `point_estimate`. `samples` is either a `NULL` or a `stanfit` object containing the samples drawn from the prior distribution of the ensemble model and `point_estimate` is either a `NULL` or a `list` object containing the optimised prior distribution of the ensemble model. If this object is `missing` then `sample_prior` generates it.
#'@param num_samples A `numeric` specifying the number of samples to be generated. The default is 1.
#'@param ... Additional arguments passed to the function \code{rstan::sampling} or  \code{rstan::optimizing}.
#'@inherit generate_sample details
#'@return An `EnsembleSample` object.
#'@references J. Durbin, S. J. Koopman (2002) A simple and efficient simulation smoother for state space time series analysis Biometrika, Volume 89, Issue 3, August 2002, Pages 603-616,
#'@references Chris M.Strickland, Ian. W.Turner, RobertDenhamb, Kerrie L.Mengersena. Efficient Bayesian estimation of multivariate state space models Computational Statistics & Data Analysis Volume 53, Issue 12, 1 October 2009, Pages 4116-4125
#'@seealso \code{\linkS4class{EnsembleFit}}, \code{\link{EnsembleSample}}, \code{\link{generate_sample}}, \code{\link{prior_ensemble_model}}
#'@export
#'@examples
#'\donttest{
#' priors <- EnsemblePrior(4)
#' prior_density <- prior_ensemble_model(priors, M = 4)
#' samples <- sample_prior(observations = list(SSB_obs, Sigma_obs),
#'              simulators = list(list(SSB_miz, Sigma_miz),
#'                                list(SSB_ewe, Sigma_ewe),
#'                                list(SSB_fs, Sigma_fs),
#'                                list(SSB_lm, Sigma_lm)),
#'              priors = priors,
#'              sam_priors = prior_density)
#' plot(samples) #Plot the prior predictive density.
#' }
sample_prior <- function(observations, simulators, priors, sam_priors, num_samples = 1, full_sample = TRUE,...){
  if(missing(sam_priors)){
    sam_priors <- prior_ensemble_model(priors, M=length(simulators),
                                       full_sample = full_sample, ...)
  }
  ens_data <- EnsembleData(observations, simulators, priors)
  fit_prior <- EnsembleFit(ens_data, sam_priors$samples, sam_priors$point_estimate)
  return(generate_sample(fit_prior,num_samples = num_samples))
}

