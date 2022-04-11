#'Generate samples of parameters from prior distribution
#'
#'Methods to generates samples of the parameters from the prior distribution of the ensemble model.
#'
#'@param priors An `EnsemblePrior` object specifying the prior distributions for the ensemble.
#'@param M A `numeric` that represents the number of simulators. The default is 1.
#'@param full_sample A `logical` that runs a full sampling of the prior density of the ensemble model if `TRUE`. If `FALSE`, returns the point estimate which maximises the prior density of the ensemble model.
#'@param ... Additional arguments passed to the function \code{rstan::sampling} or  \code{rstan::optimizing}.
#'@return A `list` containing two items named `samples` and `point_estimate`. If `full_sample==TRUE`, `samples` is a `stanfit` and `point_estimate` is a `NULL` object, else `samples` is a `NULL` and `point_estimate` is a `list` object.
#'@references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'@seealso \code{\linkS4class{EnsembleFit}}
#'@export
#'@examples
#'\donttest{
#' num_species <- 4
#' priors <- EnsemblePrior(
#'  d = num_species,
#'  ind_st_params = list("lkj",  list(3, 2), 3),
#'  ind_lt_params = list(
#'    "beta",
#'    list(c(10,4,8, 7),c(2,3,1, 4)),
#'    list(matrix(5, num_species, num_species),
#'         matrix(0.5, num_species, num_species))
#'  ),
#'  sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'  sha_lt_params = 5,
#'  truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#' )
#' prior_density <- prior_ensemble_model(priors, M = 3)
#' }
prior_ensemble_model <- function(priors,M=1,
                                    full_sample = TRUE, ...){
  stan_input <- priors@priors_stan_input
  stan_input$M <- M
  stan_input$N <- priors@d

  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(stanmodels$ensemble_prior, data=stan_input,...)
  }else{
    point_estimate <- rstan::optimizing(stanmodels$ensemble_prior, data=stan_input,as_vector=FALSE, ...)
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
#'@references J. Durbin, S. J. Koopman (2002) A simple and efficient simulation smoother for state space time series analysis Biometrika, Volume 89, Issue 3, August 2002, Pages 603–616,
#'@references Chris M.Strickland, Ian. W.Turner, RobertDenhamb, Kerrie L.Mengersena. Efficient Bayesian estimation of multivariate state space models Computational Statistics & Data Analysis Volume 53, Issue 12, 1 October 2009, Pages 4116-4125
#'@seealso \code{\linkS4class{EnsembleFit}}, \code{\link{EnsembleSample}}, \code{\link{generate_sample}}, \code{\link{prior_ensemble_model}}
#'@export
#'@examples
#'\donttest{
#' num_species <- 4
#' priors <- EnsemblePrior(
#'  d = num_species,
#'  ind_st_params = list("lkj",  list(3, 2), 3),
#'  ind_lt_params = list(
#'    "beta",
#'    list(c(10,4,8, 7),c(2,3,1, 4)),
#'    list(matrix(5, num_species, num_species),
#'         matrix(0.5, num_species, num_species))
#'  ),
#'  sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'  sha_lt_params = 5,
#'  truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#' )
#' prior_density <- prior_ensemble_model(priors, M = 3)
#' samples <- sample_prior(observations = list(SSB_obs, Sigma_obs),
#'              simulators = list(list(SSB_miz, Sigma_miz),
#'                                list(SSB_ewe, Sigma_ewe),
#'                                list(SSB_fs, Sigma_fs)),
#'              sam_priors = prior_density)
#' plot(samples) #Plot the prior predictive densty.
#' }
sample_prior <- function(observations, simulators, priors, sam_priors, num_samples = 1, full_sample = TRUE,...){
  if(missing(sam_priors)){
    sam_priors <- prior_ensemble_model(priors,
                         full_sample = full_sample,M=length(simulators), ...)
  }
  ens_data <- EnsembleData(observations, simulators, priors)
  fit_prior <- EnsembleFit(ens_data, sam_priors$samples, sam_priors$point_estimate)
  return(generate_sample(fit_prior,num_samples = num_samples))
}

