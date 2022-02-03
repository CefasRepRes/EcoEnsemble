#'Return the compiled ensemble model Stan object.
#'
#'Gets the unfit, compiled Stan object encoding the ensemble model. This allows for
#'manual fitting of the ensemble model directly using.
#'@return The `stanmodel` object encoding the ensemble model.
#'@export
#'@examples
#'mod <- get_mcmc_ensemble_model()
#'ensemble_data <- EnsembleData(observations = list(SSB_obs, Sigma_obs),
#'                              simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "mizer")),
#'                               priors = priors)
#'out <- rstan::sampling(mod, ensemble_data@@stan_input, chains=(parallel::detectCores()-1))

get_mcmc_ensemble_model <- function(){
  return(stanmodels$ensemble_model)
}

#'Fits the ensemble model
#'
#'`fit_ensemble_model` runs an MCMC of the ensemble model. This process can take a long time depending on the size of the datasets.
#'
#'@inheritParams EnsembleData
#'@param full_sample A `logical` that runs a full sampling of the ensemble model if `TRUE` and optimises to the MLE if `FALSE`. The default is TRUE.
#'@param ... Additional parameters passed to the function \code{rstan::sampling}.
#'@inherit EnsembleData details
#'@return An `EnsembleFit` object containing the data passed through to  and the fitted `stanfit` object.
#'@references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'@seealso \code{\linkS4class{EnsembleFit}}, \code{\link{EnsembleSample}}
#'@export
#'@examples
#' N_species <- 4
#' priors <- define_priors(ind_st_var_params = list(25, 0.25),
#'                         ind_st_cor_form = "lkj",
#'                         ind_st_cor_params = 30,
#'                         ind_lt_var_params = list(rep(25,N_species),rep(0.25,N_species)),
#'                         ind_lt_cor_form = "beta",
#'                         ind_lt_cor_params = list(matrix(40,N_species, N_species), matrix(40, N_species, N_species)),
#'
#'                         sha_st_var_exp = 3,
#'                         sha_st_cor_form = "lkj",
#'                         sha_st_cor_params = 30,
#'                         sha_lt_sd = rep(4,N_species))
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
fit_ensemble_model <- function(observations, simulators, priors,
                               full_sample = TRUE, ...){

  ens_data <- EnsembleData(observations, simulators, priors)
  stan_input <- ens_data@stan_input

  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(stanmodels$ensemble_model, data=stan_input,...)
  }else{
    point_estimate <- rstan::optimizing(stanmodels$ensemble_model, data=stan_input,
                                        as_vector=FALSE, ...)
  }

  return(EnsembleFit(ens_data, samples, point_estimate))
}
