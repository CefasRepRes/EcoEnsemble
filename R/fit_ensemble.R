#'Return the compiled ensemble model Stan object.
#'
#'Gets the unfit, compiled `stanmodel` object encoding the ensemble model. This allows for
#'manual fitting of the ensemble model directly using.
#'@return The `stanmodel` object encoding the ensemble model.
#'@export
#'@examples
#'mod <- get_mcmc_ensemble_model()
#'
#'priors <- EnsemblePrior(4)
#'ensemble_data <- EnsembleData(observations = list(SSB_obs, Sigma_obs),
#'                              simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "mizer")),
#'                               priors = priors)
#'\donttest{
#'out <- rstan::sampling(mod, ensemble_data@@stan_input, chains = 1)
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
#'@param control If creating a full sample, this is a named `list` of paramaters to control Stan's sampling behaviour. See the documentation of the `stan()` function in the `rstan` package for details. The default value is `list(adapt_delta = 0.95)`. If optimizing, this value is ignored.
#'@param ... Additional arguments passed to the function \code{rstan::sampling} or  \code{rstan::optimizing}.
#'@inherit EnsembleData details
#'@return An `EnsembleFit` object.
#'@references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'@seealso \code{\linkS4class{EnsembleFit}}, \code{\link{EnsembleSample}}
#'@export
#'@examples
#'\donttest{
#'
#' fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                  list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                  list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                  list(SSB_miz, Sigma_miz, "Mizer")),
#'                priors = EnsemblePrior(4,
#'                ind_st_params = IndSTPrior(parametrisation_form = "lkj",
#'                var_params= list(1,1), cor_params = 10, AR_params = c(2, 2))),
#'                full_sample = FALSE) #Only optimise in this case
#'}
fit_ensemble_model <- function(observations, simulators, priors,
                               full_sample = TRUE, control = list(adapt_delta = 0.95), ...){

  ens_data <- EnsembleData(observations, simulators, priors)
  stan_input <- ens_data@stan_input

  #Using hierarchical priors uses a different model. This speeds up the sampling enormously
  mod <- stanmodels$ensemble_model
  if(stan_input$form_prior_ind_st == 3){
    mod <- stanmodels$ensemble_model_hierarchical
    if(!full_sample){
      stop("It is possible to generate a point estimate for the prior if the individual short-term discrepancy prior follows a hierarchical parameterisation. Please generate a full sample using 'full_sample=TRUE'.")
    }
  }


  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(mod, data=stan_input, control = control, ...)
  }else{
    point_estimate <- rstan::optimizing(mod, data=stan_input,as_vector=FALSE, ...)
  }

  return(EnsembleFit(ens_data, samples, point_estimate))
}
