#Returns the compiled ensemble model Stan object in the case with drivers.


fit_ensemble_model_dri <- function(observations, simulators, priors,
                               full_sample = TRUE, control = list(adapt_delta = 0.95), MMod, ...){


  if (!missing(MMod)) {
    #ens_data <- EnsembleDataDriver(observations, simulators, priors, MMod)
    ens_data <- EnsembleData(observations, simulators, priors, drivers = TRUE, MMod)
  }
  else {
    #ens_data <- EnsembleDataDriver(observations, simulators, priors)
    ens_data <- EnsembleData(observations, simulators, priors, drivers = TRUE)
  }
  stan_input <- ens_data@stan_input

  #Using hierarchical priors uses a different model. This speeds up the sampling enormously
  mod <- stanmodels$ensemble_model_withdrivers
  if(stan_input$form_prior_ind_st == 3 || stan_input$form_prior_ind_st == 4){
    mod <- stanmodels$ensemble_model_hierarchical_withdrivers
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
