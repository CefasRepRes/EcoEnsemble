#'Generate samples from a fitted ensemble model to get MCMC effective sample size diagnostics.
#'
#'Methods to generates samples of the latent variables from a fitted ensemble model in the same order as the MCMC and to calculate diagnostics effective sample size for the MCMC.
#'@param fit An `EnsembleFit` object.
#'@param ex.fit A named `numeric`.
#'@param transformed_data A `list` of transformed input data.
#'@param time A `numeric` specifying the time for which the ensemble model was run.
#'@param arr A `named numeric`
#'@param stan_name A `character` object
#'@param only_voi A `logical` indicating whether only the effective sample size of the variables of interest should be returned. Default value is `TRUE`.
#'@param sammy An `array` of size the number of parameters \eqn{\times} `num_samples` \eqn{\times} `nchains`.
#'
#'@details  Effective sample size is calculated from these samples, only for the quantity of interest if desired.
#'
#'@return
#'`generate_sample_array` gives an array of size  `time * (MM + M + 2)*N` \eqn{\times} `num_samples` \eqn{\times} `nchains` where `M` is the number of simulators, `MM` is the number of drivers (usually 1), `N` is the number of variables of interest `num_samples` is the number of samples from the ensemble model for each chain and `nchains` is the number of chains in the `EnsembleFit` object.
#'
#'`get_parameters_array` gives a `list` of ensemble parameters from the requested sample.
#'
#'`get_mle_array` gives an array of size `time * (MM + M + 2)*N` \eqn{\times} `num_samples` \eqn{\times} `nchains`.
#'
#'`gen_sample_array` gives a `list` of two arrays of dimensions `time * (MM + M + 2)*N` \eqn{\times} `num_samples` \eqn{\times} `nchains`
#'
#'`get_ESS_diag` and `calc_ess` gives a list of length two containing `ESS_bulk` and `ESS_tail` which has the bulk effective sample size and the tail effective sample size. Each list is of dimensions `time` \eqn{\times} `(MM + M + 2)*N`.
#'
#' gives a list of length two containing `ESS_bulk` and `ESS_tail` which has the bulk effective sample size and the tail effective sample size.
#'@importFrom posterior ess_bulk
#'@importFrom posterior ess_tail
#'
#'@rdname get_stan_outputs_array
#'@seealso See [generate_sample()] for information about the sampling and \code{\link[posterior]{ess_bulk}} and \code{\link[posterior]{ess_tail}} for information about the effective sample size calculations.
#'@export
#'@examples
#'\donttest{
#' fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                  list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                  list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                  list(SSB_miz, Sigma_miz, "Mizer")),
#'                priors = EnsemblePrior(4,
#'                ind_st_params = IndSTPrior(parametrisation_form = "lkj",
#'                var_params= list(1,1), cor_params = 10, AR_params = c(2, 2))))
#' get_ESS_diag(fit)
#'}
generate_sample_array <- function(fit) {
  stan_input <- fit@ensemble_data@stan_input

  # If we have samples, then a full MCMC was run
  full_sample <- !is.null(fit@samples)

  if ("MM" %in% names(stan_input)) {
    transformed_data <- get_transformed_data_dri(fit)
  } else {
    transformed_data <- get_transformed_data(fit)
  }

  if (!full_sample) {
    stop("`fit@samples` is NULL. A full MCMC fit is required.", call. = FALSE)
  }

  ex.fit <- rstan::extract(fit@samples, permuted = FALSE)


  mle <- apply(
    ex.fit,
    c(1, 2),
    get_mle_array,
    transformed_data = transformed_data,
    time = stan_input$time
  )

  sammy <- apply(
    ex.fit,
    c(1, 2),
    gen_sample_array,
    transformed_data = transformed_data,
    time = stan_input$time
  )

  hat <- array(unlist(lapply(sammy, "[[", "sam_x_hat")), dim = dim(mle))
  x <- array(unlist(lapply(sammy, "[[", "sam_x")), dim = dim(mle))

  sample_ret <- mle - hat + x

  return(sample_ret)
}

#'@rdname get_stan_outputs_array
#'@export
get_mle_array <- function(ex.fit, transformed_data, time) {
  params <- get_parameters_array(ex.fit)

  ret <- KalmanFilter_back(
    params$AR_params,
    params$lt_discrepancies,
    transformed_data$all_eigenvalues_cov,
    params$SIGMA,
    transformed_data$bigM,
    params$SIGMA_init,
    params$x_hat,
    time,
    transformed_data$new_data,
    transformed_data$observation_available
  )

  return(ret)
}

#'@rdname get_stan_outputs_array
#'@export
get_parameters_array <- function(ex.fit) {
  ret <- list()

  tmp <- ex.fit[get_param_idx(ex.fit, "x_hat")]
  dimmy <- length(tmp)

  ret$x_hat <- tmp
  ret$SIGMA_init <- matrix(
    ex.fit[get_param_idx(ex.fit, "SIGMA_init")],
    dimmy,
    dimmy
  )
  ret$AR_params <- ex.fit[get_param_idx(ex.fit, "AR_params")]
  ret$SIGMA <- matrix(
    ex.fit[get_param_idx(ex.fit, "SIGMA")],
    dimmy,
    dimmy
  )
  ret$lt_discrepancies <- ex.fit[get_param_idx(ex.fit, "lt_discrepancies")]

  return(ret)
}

#'@rdname get_stan_outputs_array
#'@export
get_param_idx <- function(arr, stan_name) {
  pnames <- names(arr)
  grep(paste0("^", stan_name, "(\\[|$)"), pnames)
}

#'@rdname get_stan_outputs_array
#'@export
gen_sample_array <- function(ex.fit, transformed_data, time) {
  params <- get_parameters_array(ex.fit)

  bigM <- transformed_data$bigM
  all_eigenvalues_cov <- transformed_data$all_eigenvalues_cov
  observation_available <- transformed_data$observation_available

  sam_y <- matrix(0, time, nrow(bigM))
  sam_x <- matrix(0, time, length(params$x_hat))

  sam_x[1, ] <- as.matrix(params$x_hat) +
    t(chol(params$SIGMA_init)) %*% rnorm(length(params$x_hat))
  sam_y[1, ] <- bigM %*% (sam_x[1, ] + params$lt_discrepancies) +
    rnorm(nrow(bigM), 0, sqrt(all_eigenvalues_cov))

  ch_sigma <- t(chol(params$SIGMA))

  for (i in 2:time) {
    sam_x[i, ] <- as.matrix(params$AR_params * sam_x[i - 1, ]) +
      ch_sigma %*% rnorm(length(params$x_hat))
    sam_y[i, ] <- bigM %*% (sam_x[i, ] + params$lt_discrepancies) +
      rnorm(nrow(bigM), 0, sqrt(all_eigenvalues_cov))
  }

  sam_x_hat <- KalmanFilter_back(
    params$AR_params,
    params$lt_discrepancies,
    all_eigenvalues_cov,
    params$SIGMA,
    bigM,
    params$SIGMA_init,
    params$x_hat,
    time,
    sam_y,
    observation_available
  )

  return(list(sam_x = sam_x, sam_x_hat = sam_x_hat))
}

#'@rdname get_stan_outputs_array
#'@export
get_ESS_diag <- function(fit,only_voi=TRUE){
  sammy <- generate_sample_array(fit)
  ####
  if (only_voi==TRUE){
    idx <- 1:(fit@ensemble_data@stan_input$time * fit@ensemble_data@stan_input$N)
    #sammy_save <- sammy
    sammy <- sammy[idx,,]
  }
  ESS_tmp <- calc_ess(sammy)
  ret <- list(
    ESS_bulk = matrix(ESS_tmp$ESS_bulk,nrow=fit@ensemble_data@stan_input$time),
    ESS_tail = matrix(ESS_tmp$ESS_tail,nrow=fit@ensemble_data@stan_input$time)
  )
  return(ret)
}

#'@rdname get_stan_outputs_array
#'@export
calc_ess <- function(sammy){
  ret <- list(
    ESS_bulk = apply(sammy,1,posterior::ess_bulk),
    ESS_tail = apply(sammy,1,posterior::ess_tail)
  )
  return(ret)
}
