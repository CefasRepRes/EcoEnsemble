#'Return the compiled ensemble model STAN object.
#'
#'Gets the unfit, compiled STAN object encoding the ensemble model. This allows for
#'manual fitting of the ensemble model directly using.
#'@return The `stanmodel` object encoding the ensemble model.
#'@export
#'@examples
#'#A minimal example fitting manually
#'mod <- get_mcmc_ensemble_model()
#'stan_data <- get_stan_ensemble_data(observations = list(SSB_obs, Sigma_obs),
#'                                    simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                                      list(SSB_fs,  Sigma_fs),
#'                                                      list(SSB_lm,  Sigma_lm),
#'                                                      list(SSB_miz, Sigma_miz)),
#'                                    priors = priors)
#'out <- rstan::sampling(mod, stan_data, chains=(parallel::detectCores()-1))

get_mcmc_ensemble_model <- function(){
  return(stanmodels$ensemble_model)
}
#'Validate and generate STAN data from inputs.
#'
#'Generates data in the correct form to be passed through to the compiled STAN model.
#'@param observations A `list` of observations and associated covariances.
#'See details.
#'@param simulators A `list` of simulator output and covariance pairs. See details.
#'@param priors A `list` of data specifying the prior distributions for the ensemble.
#'See the `elicit_priors` function for the required form.
#'@details Observation and covariance pairs should be passed through as a `list`,
#'with observations / model outputs as the first element, and the covariance matrix
#'as the second. The observations / model outputs should be a matrix / data frame with each
#'column giving a different variable / species and each row a different year. Rows should be named with
#'the years and columns should be named the appropriate variables. It is fine to have missing
#'years or species, however for any model output there should be at least one observation in the observation
#'data frame.
#'@return An `EnsembleData` object containing the data to be passed through to STAN.
#'@export
get_stan_ensemble_data <- function(observations, simulators, priors){

  validate_data(observations, simulators, priors)

  stan_data_sims <- generate_simulator_stan_data(observations, simulators)
  stan_data_priors <- priors[- which(names(priors) == "prior_correlations" )]
  stan_data_priors <- append(stan_data_priors,
                             generate_correlation_priors_stan_data(priors$prior_correlations))

  ret <- append(stan_data_sims, stan_data_priors)
  class(ret) <- "EnsembleData"

  return(ret)
}

#'Fits the ensemble model
#'
#'`fit_ensemble_model` runs an MCMC of the ensemble model through STAN. This process can take a long time depending on
#'the size of the datasets.
#'
#'@inheritParams get_stan_ensemble_data
#'@param full_sample A `logical` that runs a full sampling of the ensemble model if `TRUE`,
#'and optimises to the MLE if `FALSE`.
#'@param ... Additional STAN parameters.
#'@inherit get_stan_ensemble_data details
#'@return An `EnsembleFit` object containing the data passed through to STAN
#'and the fitted `stanfit` object.
#'
#'@export
#'@examples
#' A basic example
#' N_species <- 4
#' priors <- elicit_priors(ind_st_var_params = list(25, 0.25),
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
#'fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                      list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                      list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                      list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors,
#'                          full_sample = FALSE, #Only optimise in this case
#'                          control = list(adapt_delta = 0.99))#Additional STAN options
fit_ensemble_model <- function(observations, simulators, priors,
                               full_sample = FALSE, ...){
  stan_input <- get_stan_ensemble_data(observations,simulators,priors)

  if(full_sample){
    ex.fit <- rstan::sampling(stanmodels$ensemble_model, data=stan_input,...)
  }else{
    ex.fit <- rstan::optimizing(stanmodels$ensemble_model, data=stan_input,
                                as_vector=FALSE, ...)

  }

  ret <- list(observations = observations, simulators = simulators, priors = priors,
              stan_input = stan_input, ex.fit = ex.fit, full_sample = full_sample)

  class(ret) <- "EnsembleFit"
  return(ret)
}
