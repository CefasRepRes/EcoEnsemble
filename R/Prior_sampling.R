prior_ensemble_model <- function(priors,M,
                                    full_sample = TRUE, ...){
  stan_input <- priors@priors_stan_input
  stan_input$M <- M

  samples <- NULL; point_estimate <- NULL
  if(full_sample){
    samples <- rstan::sampling(stanmodels$ensemble_prior, data=stan_input,...)
  }else{
    point_estimate <- rstan::optimizing(stanmodels$ensemble_prior, data=stan_input,as_vector=FALSE, ...)
  }
  return(list(samples=samples, point_estimate=point_estimate))
}

sample_prior <- function(observations, simulators, priors,sam_priors,num_samples = 1,...){
  if(is.missing(sam_priors)){
    sam_priors <- prior_ensemble_model(priors,
                         full_sample = TRUE,M=length(simulators), ...)
  }
  ens_data <- EnsembleData(observations, simulators, priors)
  fit_prior <- EnsembleFit(ens_data, sam_priors$samples, sam_priors$point_estimate)
  return(generate_sample(fit_prior,num_samples = num_samples))
}

