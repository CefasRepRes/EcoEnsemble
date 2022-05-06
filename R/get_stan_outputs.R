#'Generate samples from a fitted ensemble model.
#'
#'Methods to generates samples of the latent variables from a fitted ensemble model.
#'@param fit An `EnsembleFit` object.
#'@param num_samples A `numeric` specifying the number of samples to be generated. The default is 1.
#'@param time A `numeric` specifying the time for which the ensemble model was run.
#'@param ex.fit A `list` containing the samples / point estimate from the `EnsembleFit` object.
#'@param x A `numeric` specifying which sample from the posterior to use. The default is 1.
#'@param transformed_data A `list` of transformed input data.
#'
#'@details The samples are created using the methods described in Strickland et. al. (2009) and Durbin and Koopman (2002).
#'
#'@return
#'`generate_sample` gives a `list` of length 2, the first element being the MLE of latent variables and the second element being a set of samples of the latent variables.
#'
#'* If `fit` is a sampling of the ensemble model parameters, then:
#'    + `mle` is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`, where \eqn{M} is the number of simulators and `num_samples` is the number of samples from the ensemble model, giving the MLE of the latent variables for each available sample from the ensemble model.
#'    + `sample` is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`, giving a sample of the latent variables for each available sample of the ensemble model.
#'
#'* If `fit` is a point estimate of the ensemble model parameters, then:
#'    + `mle` is a `time`\eqn{\times (3M + 2) \times} 1 `array` giving the MLE of the latent variables for the point estimate of the ensemble model.
#'    + `sample` is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`, giving `num_samples` samples of the latent variables for the single point estimate of the ensemble model.
#'
#'`get_transformed_data` gives a `list` of transformed input data.
#'
#'`get_parameters` gives a `list` of ensemble parameters from the requested sample.
#'
#'`get_mle` If `fit` is a sampling of the ensemble model parameters, then this is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`. If `fit` is a point estimate of the ensemble model parameters, then this is a `time`\eqn{\times (3M + 2) \times} 1 `array` giving the MLE of the latent variables for the point estimate of the ensemble model.
#'
#'`gen_sample` If `fit` is a sampling of the ensemble model parameters, then this is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`, giving a sample of the latent variables for each available sample of the ensemble model. If `fit` is a point estimate of the ensemble model parameters, then this is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`.
#'
#'
#'@rdname get_stan_outputs
#' @references J. Durbin, S. J. Koopman (2002) A simple and efficient simulation smoother for state space time series analysis Biometrika, Volume 89, Issue 3, August 2002, Pages 603–616,
#' @references Chris M.Strickland, Ian. W.Turner, RobertDenhamb, Kerrie L.Mengersena. Efficient Bayesian estimation of multivariate state space models Computational Statistics & Data Analysis Volume 53, Issue 12, 1 October 2009, Pages 4116-4125
#'@export
#'@examples
#'num_species <- 4
#' priors <- EnsemblePrior(
#'     d = num_species,
#'     ind_st_params = list("lkj",  list(3, 2), 3),
#'     ind_lt_params = list(
#'        "beta",
#'         list(c(10,4,8, 7),c(2,3,1, 4)),
#'         list(matrix(5, num_species, num_species),
#'              matrix(0.5, num_species, num_species))
#'     ),
#'     sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'     sha_lt_params = 5,
#'     truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#' )
#'\donttest{
#'fit <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors)
#'sample <- generate_sample(fit, num_samples = 2000)
#'
#'# A quicker way to get the MLE for the first sample from the ensemble
#'transf_data <- get_transformed_data(fit)
#'ex.fit <- rstan::extract(fit@@samples)
#'mle_sample <- get_mle(1, ex.fit = ex.fit, transformed_data = transformed_data,
#'                      time = fit@@ensemble_data@@stan_input$time, simplify = F)
#'}
generate_sample <- function(fit, num_samples = 1)
{

  stan_input <- fit@ensemble_data@stan_input

  # If we have samples, then a full MCMC was run
  full_sample <- !is.null(fit@samples)

  transformed_data <- get_transformed_data(fit)

  rstan::expose_stan_functions(stanmodels$KF_back)

  if (full_sample){
    ex.fit <- rstan::extract(fit@samples)
    num_samples <- nrow(ex.fit$x_hat)
    mle <- sapply(1:num_samples, get_mle, ex.fit = ex.fit, transformed_data = transformed_data,
                      time = stan_input$time, simplify = F)
  }else {
    ex.fit <- fit@point_estimate$par
    mle <- get_mle(x=1, ex.fit, transformed_data = transformed_data,
                       time = stan_input$time)
  }

  sammy <- sapply(1:num_samples, gen_sample, ex.fit = ex.fit,
                  transformed_data = transformed_data,
                  time = stan_input$time, simplify = F)


  if (full_sample){
    tmp <- mapply(function(x,y){
      y - x$sam_x_hat + x$sam_x
    }, y = mle,
    x = sammy)
    sample_ret <-array(as.numeric(unlist(tmp)),dim=c(nrow(mle[[1]]),ncol(mle[[1]]),num_samples))
    mle <-array(as.numeric(unlist(mle)),dim=c(nrow(mle[[1]]),ncol(mle[[1]]),num_samples))
  }else{
    tmp<-lapply(sammy,function(x,y){
      y - x$sam_x_hat + x$sam_x
    }, y = mle)
    sample_ret <-array(as.numeric(unlist(tmp)),dim=c(nrow(mle), ncol(mle), num_samples))
  }

  return(EnsembleSample(fit, mle, sample_ret))
}


# Fudge to get past "no visible binding for global variable" in R-CMD check
utils::globalVariables(c("KalmanFilter_back"))

#'@rdname get_stan_outputs
#'@export
get_transformed_data <- function(fit){
  stan_input <- fit@ensemble_data@stan_input

  N <- stan_input$N
  M <- stan_input$M
  model_num_species <- stan_input$model_num_species
  obs_covariances <- stan_input$obs_covariances
  Ms <- stan_input$Ms
  model_covariances <- stan_input$model_covariances
  observation_times <- stan_input$observation_times
  time <- stan_input$time
  observations <- stan_input$observations
  model_outputs <- stan_input$model_outputs

  bigM <- matrix(0,N + sum(model_num_species) , (M+2) * N )
  all_eigenvalues_cov <- rep(0,N + sum(model_num_species) )
  all_eigenvectors_cov = matrix(0,N + sum(model_num_species) , N + sum(model_num_species) )
  new_data <- observation_available <- matrix(0,time , N + sum(model_num_species))
  bigM[1:N,1:N] = diag(N)

  tmp_eigen <- eigen(obs_covariances,symmetric=T)
  all_eigenvalues_cov[1:N] <- tmp_eigen$values
  all_eigenvectors_cov[1:N,1:N] <- tmp_eigen$vectors
  if (M > 0){
    bigM[(N + 1):(N + model_num_species[1] ),1:N] <- Ms[1:model_num_species[1],]
    bigM[(N + 1):(N + model_num_species[1] ),(N+1):(2*N)] <-  Ms[1:model_num_species[1],]
    bigM[(N + 1):(N + model_num_species[1] ),(2*N + 1):(3 * N)] <- Ms[1:model_num_species[1],]

    tmp_eigen <- eigen(matrix(model_covariances[1:(model_num_species[1]^2)],model_num_species[1],model_num_species[1]),symmetric=T)
    all_eigenvalues_cov[(N + 1):(N + model_num_species[1] )] <- tmp_eigen$values
    all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] = tmp_eigen$vectors;

    if(M > 1){
      for(i in 2:M){
        start_index <-  N + sum(model_num_species[1:(i-1)]) + 1
        end_index <- N + sum(model_num_species[1:i])
        Ms_transformation <- Ms[(sum(model_num_species[1:(i-1)]) + 1):sum(model_num_species[1:i]),]
        bigM[start_index:end_index,1:N] <- Ms_transformation
        bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation
        bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] <- Ms_transformation

        tmp_eigen <- eigen(matrix(model_covariances[(cumsum(model_num_species^2)[i-1] + 1):(cumsum(model_num_species^2)[i])],model_num_species[i],model_num_species[i]),symmetric=T)
        all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
        all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors
      }
    }
  }

  for (i in 1:time){
    observation_available[i,1:N] = rep(observation_times[i,1],N)
    if (M > 0){
      observation_available[i,(N+1):(N + model_num_species[1])] <- rep(observation_times[i,2],model_num_species[1])
      if(M > 1){
        for (j in 2:M){
          observation_available[i,(N + sum(model_num_species[1:(j-1)]) + 1):(N + sum(model_num_species[1:j]))] <- rep(observation_times[i,j + 1],model_num_species[j])
        }
      }
    }

    new_data[i,] <- t(all_eigenvectors_cov) %*% matrix(c(observations[i,],model_outputs[i,]),ncol=1)
  }

  bigM <- t(all_eigenvectors_cov) %*% bigM

  return(
    list(bigM=bigM,
         all_eigenvalues_cov=all_eigenvalues_cov,
         new_data=new_data,
         observation_available=observation_available)
  )
}

#'@rdname get_stan_outputs
#'@export
get_parameters <- function(ex.fit, x = 1){
  ret <- list()
  if (is(ex.fit$x_hat,"matrix")){ ## It's a sample
    ret$x_hat <- ex.fit$x_hat[x,]
    ret$SIGMA_init <- ex.fit$SIGMA_init[x,,]
    ret$AR_params <- ex.fit$AR_params[x,]
    ret$new_x <- ex.fit$x_hat[x,]
    ret$SIGMA <- ex.fit$SIGMA[x,,]
    ret$lt_discrepancies <- ex.fit$lt_discrepancies[x,]
  }else{ ## It's optimised
    ret$x_hat <- ex.fit$x_hat
    ret$SIGMA_init <- ex.fit$SIGMA_init
    ret$AR_params <- ex.fit$AR_params
    ret$new_x <- ex.fit$x_hat
    ret$SIGMA <- ex.fit$SIGMA
    ret$lt_discrepancies <- ex.fit$lt_discrepancies
  }
  return(ret)
}


#'@rdname get_stan_outputs
#'@export
get_mle <- function(x=1, ex.fit, transformed_data, time) ## get the MLE from the fit
{

  if(!exists("KalmanFilter_back"))
   rstan::expose_stan_functions(stanmodels$KF_back)

  params <- get_parameters(ex.fit,x)

  ret <- KalmanFilter_back(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                           params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                           transformed_data$new_data, transformed_data$observation_available)
  return(ret)
}


#'@rdname get_stan_outputs
#'@export
gen_sample <- function(x=1, ex.fit, transformed_data, time)
{
  if(!exists("KalmanFilter_back"))
    rstan::expose_stan_functions(stanmodels$KF_back)

  params <- get_parameters(ex.fit, x)

  bigM = transformed_data$bigM
  all_eigenvalues_cov = transformed_data$all_eigenvalues_cov
  observation_available = transformed_data$observation_available

  sam_y <- matrix(0,time,nrow(bigM))
  sam_x <- matrix(0,time,length(params$x_hat))

  sam_x[1,] <- as.matrix(params$x_hat) + chol(params$SIGMA_init)%*%rnorm(length(params$x_hat))
  sam_y[1,] <- bigM %*% (sam_x[1,] + params$lt_discrepancies) + rnorm(nrow(bigM),0,sqrt(all_eigenvalues_cov))
  ch_sigma <- chol(params$SIGMA)
  for (i in 2:time){
    sam_x[i,] <- as.matrix(params$AR_params * sam_x[i-1,]) + ch_sigma%*%rnorm(length(params$x_hat))
    sam_y[i,] <- bigM %*% (sam_x[i,] + params$lt_discrepancies) + rnorm(nrow(bigM),0,sqrt(all_eigenvalues_cov))
  }
  sam_x_hat <- KalmanFilter_back(params$AR_params, params$lt_discrepancies,
                                 all_eigenvalues_cov, params$SIGMA,
                                 bigM,
                                 params$SIGMA_init,
                                 params$x_hat, time, sam_y,
                                 observation_available)
  return(list(sam_x=sam_x,sam_x_hat=sam_x_hat))
}
