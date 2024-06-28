#'Return the compiled ensemble model Stan object.
#'
#'Gets the unfit, compiled `stanmodel` object encoding the ensemble model. This allows for
#'manual fitting of the ensemble model directly using \code{rstan::sampling}.
#'
#'@param priors An `EnsemblePrior` object specifying the prior distributions for the ensemble for which the compiled `stanmodel` object will be obtained.
#'@param likelihood A `logical` that returns the compiled `stanmodel` object including the likelihood (the Kalman filter) for given priors if `TRUE`. If `FALSE` returns the compiled `stanmodel` object without the likelihood for sampling from the prior.
#'@param drivers A `logical` indicating whether drivers have been used in combination with simulators. Default value is FALSE.
#'@return The `stanmodel` object encoding the ensemble model.
#'@export
#'@examples
#'priors <- EnsemblePrior(4)
#'mod <- get_mcmc_ensemble_model(priors)
#'
#'ensemble_data <- EnsembleData(observations = list(SSB_obs, Sigma_obs),
#'                              simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "mizer")),
#'                               priors = priors)
#'\donttest{
#'out <- rstan::sampling(mod, ensemble_data@@stan_input, chains = 1)
#'}

get_mcmc_ensemble_model <- function(priors, likelihood = TRUE, drivers = FALSE){

  st_pf <- priors@ind_st_params@parametrisation_form

  if (st_pf == "hierarchical") {
    if (likelihood) {
      if (!drivers){
        return(stanmodels$ensemble_model_hierarchical)
      } else {
        return(stanmodels$ensemble_model_hierarchical_withdrivers)
      }
    } else {
      if (!drivers) {
        return(stanmodels$ensemble_prior_hierarchical)
      } else {
        return(stanmodels$ensemble_prior_hierarchical_withdrivers)
      }
    }
  } else {
    if (likelihood) {
      if (!drivers) {
        return(stanmodels$ensemble_model)
      } else {
        return(stanmodels$ensemble_model_withdrivers)
      }
    } else {
      if (!drivers) {
        return(stanmodels$ensemble_prior)
      } else {
        return(stanmodels$ensemble_prior_withdrivers)
      }
    }
  }
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
                               full_sample = TRUE, control = list(adapt_delta = 0.95), drivers = FALSE, MMod, ...){

  if (drivers == TRUE) {
    return(fit_ensemble_model_dri(observations, simulators, priors,
                                  full_sample = full_sample, control = list(adapt_delta = 0.95), MMod, ...))
  }
  else {
    ens_data <- EnsembleData(observations, simulators, priors)
    stan_input <- ens_data@stan_input

    #Using hierarchical priors uses a different model. This speeds up the sampling enormously
    mod <- stanmodels$ensemble_model
    if(stan_input$form_prior_ind_st == 3 || stan_input$form_prior_ind_st == 4){
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
}
