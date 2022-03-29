#'Return the compiled ensemble model Stan object.
#'
#'Gets the unfit, compiled Stan object encoding the ensemble model. This allows for
#'manual fitting of the ensemble model directly using.
#'@return The `stanmodel` object encoding the ensemble model.
#'@export
#'@examples
#'mod <- get_mcmc_ensemble_model()
#'
#'num_species <- 4
#'priors <- EnsemblePrior(
#'  d = num_species,
#'  ind_st_params = list("lkj",  list(3, 2), 3),
#'  ind_lt_params = list(
#'     "beta",
#'     list(c(10,4,8, 7),c(2,3,1, 4)),
#'     list(matrix(5, num_species, num_species),
#'          matrix(0.5, num_species, num_species))
#'  ),
#'  sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'  sha_lt_params = 5,
#'  truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#')
#'ensemble_data <- EnsembleData(observations = list(SSB_obs, Sigma_obs),
#'                              simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "mizer")),
#'                               priors = priors)
#'\dontrun{
#'out <- rstan::sampling(mod, ensemble_data@@stan_input, chains=(parallel::detectCores()-1))
#'}

get_mcmc_ensemble_model <- function(){
  return(stanmodels$ensemble_model)
}

#'Fits the ensemble model
#'
#'`fit_ensemble_model` runs an MCMC of the ensemble model. This process can take a long time depending on the size of the datasets.
#'
#'@inheritParams EnsembleData
#'@param full_sample A `logical` that runs a full sampling of the posterior density of the ensemble model if `TRUE`. If `FALSE`, returns the point estimate which maximises the posterior density of the ensemble model.
#'@param ... Additional arguments passed to the function \code{rstan::sampling} or  \code{rstan::optimizing}.
#'@inherit EnsembleData details
#'@return An `EnsembleFit` object.
#'@references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'@seealso \code{\linkS4class{EnsembleFit}}, \code{\link{EnsembleSample}}
#'@export
#'@examples
#'\dontrun{
#' num_species <- 4
#' priors <- EnsemblePrior(
#'     d = num_species,
#'     ind_st_params = list("lkj",  list(3, 2), 3),
#'     ind_lt_params = list(
#'        "beta",
#'        list(c(10,4,8, 7),c(2,3,1, 4)),
#'        list(matrix(5, num_species, num_species),
#'             matrix(0.5, num_species, num_species))
#'     ),
#'     sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'     sha_lt_params = 5,
#'     truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#' )
#' fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                      list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                      list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                      list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors,
#'                          full_sample = FALSE) #Only optimise in this case
#' Run a full sample
#' fit1 <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                      list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                      list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                      list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors,
#'                          control=list(adapt_delta = 0.99)) # Additional Stan parameters.
#'}
fit_ensemble_model <- function(observations, simulators, priors,
                               full_sample = TRUE, ...){

  ens_data <- EnsembleData(observations, simulators, priors)
  stan_input <- ens_data@stan_input

  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(stanmodels$ensemble_model, data=stan_input,...)
  }else{
    point_estimate <- rstan::optimizing(stanmodels$ensemble_model, data=stan_input,as_vector=FALSE, ...)
  }

  return(EnsembleFit(ens_data, samples, point_estimate))
}
