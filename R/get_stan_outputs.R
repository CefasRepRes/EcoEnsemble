#'Generate samples from a fitted ensemble model.
#'
#'Methods to generates samples of the latent variables in the ensemble model by sampling from
#'the Kalman filter. It is strongly recommended to use the `generate_sample` wrapper function,
#'but subroutiones are provided to allow finer control if desired.
#'@param fit An `EnsembleFit` fitted ensemble object, an output from the `fit_ensemble_model` function.
#'@param num_samples The number of samples to generate from the Kalman filter if the `EnsembleFit` object is not a full sampling
#'of the posterior of the ensemble model. The default is 1. If the `EnsembleFit` object is a full sampling of the posterior, then the
#'number of samples will be the same as the number of samples in the ensemble MCMC.
#'@param time A non-negative integer specifying the time for which the ensemble model was run.
#'@param ex.fit The extracted samples / point estimate from the `EnsembleFit` object. For samples, this should be the output
#'of `rstan::extract()` on the `stanfit` object and for point estimate outputs this should be the `par` slot.
#'@param x The sample of parameters from the ensemble model to use (if the `EnsembleFit` object is a full sampling, otherwise
#'uses the MLE). The default is 1.
#'
#'@details
#'
#'@return
#'`generate_sample` gives a `list` of length 2, with the first element being the MLE of the Kalman filter and
#'the second element being samples.
#'
#'* If `fit` is a full MCMC sampling of the ensemble model, then:
#'    + `mle` is a `time`\eqn{\times (3M + 2) \times N_{ens}} `array` (where \eqn{M} is the number of simulators
#' and \eqn{N_{ens}} is the number of samples from the ensemble model), giving the MLE of the Kalman filter for
#' each available sample of the ensemble model.
#'    + `sample` is a `time`\eqn{\times (3M + 2) \times N_{ens}} `array`, giving a sample of the Kalman filter for
#' each available sample of the ensemble model.
#'
#'* If `fit` is a point estimate of the ensemble model, then:
#'    + `mle` is a `time`\eqn{\times (3M + 2) \times} 1 `array` giving the MLE of the Kalman filter for the point estimate
#' of the ensemble model.
#'    + `sample` is a `time`\eqn{\times (3M + 2) \times} `num_samples` `array`, giving `num_samples` samples of the Kalman filter for
#' the single point estimate of the ensemble model.
#'
#'`get_transformed_data` gives a `list` of transformed input data in discrepancy space.
#'
#'`get_parameters` gives a `list` of ensemble parameters from the requested sample.
#'
#'`get_mle` gives the MLE of the Kalman filter for the the requested sample of ensemble parameter values.
#'
#'`gen_sample` gives a sample from the Kalman filter for the the requested sample of ensemble parameter values.
#'
#'
#'@rdname get_stan_outputs
#'@export
#'@examples
#'# Basic usage of the generate sample function
#'fit_sample <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors)
#'sample <- generate_sample(fit_sample, num_samples = 2000)
#'
#'# A quicker way to get the MLE for the first sample from the ensemble
#'transf_data <- get_transformed_data(fit_sample)
#'ex.fit <- rstan::extract(fit$ex.fit)
#'mle_sample <- get_mle(1, ex.fit = ex.fit, transformed_data = transformed_data,
#'                      time = fit$stan_input$time, simplify = F)
#'
generate_sample <- function(fit, num_samples = 1)
{

  stan_input <- fit$stan_input

  transformed_data <- get_transformed_data(fit)

  ret <- list(mle=c(),sample=c())

  rstan::expose_stan_functions(stanmodels$KF_back)

  if (fit$full_sample){
    ex.fit <- rstan::extract(fit$ex.fit)
    num_samples <- nrow(ex.fit$x_hat)
    #ret$mle <- sapply(1:num_samples, get_mle, ex.fit = ex.fit, bigM = bigM,
    #                  all_eigenvalues_cov = all_eigenvalues_cov,
    #                  observation_available = observation_available,
    #                  time = stan_input$time, new_data = new_data, simplify = F)
    ret$mle <- sapply(1:num_samples, get_mle, ex.fit = ex.fit, transformed_data = transformed_data,
                      time = stan_input$time, simplify = F)
  }else {
    ex.fit <- fit$ex.fit$par
    #ret$mle <- get_mle(x=1, ex.fit, bigM, all_eigenvalues_cov,
    #                   observation_available, time = stan_input$time, new_data)
    ret$mle <- get_mle(x=1, ex.fit, transformed_data = transformed_data,
                       time = stan_input$time)
  }

  #sammy <- sapply(1:num_samples, gen_sample, ex.fit = ex.fit, bigM=bigM,
  #                all_eigenvalues_cov = all_eigenvalues_cov,
  #                observation_available = observation_available,
  #                time = stan_input$time, simplify = F)
  sammy <- sapply(1:num_samples, gen_sample, ex.fit = ex.fit,
                  transformed_data = transformed_data,
                  time = stan_input$time, simplify = F)


  if (fit$full_sample){
    tmp <- mapply(function(x,y){
      y - x$sam_x_hat + x$sam_x
    }, y = ret$mle,
    x = sammy)
    #TODO: CHange the dimnames to be something helpful for the end user....
    #ret$sample <-array(as.numeric(unlist(tmp)),dim=c(nrow(ret$mle[[1]]),ncol(ret$mle[[1]]),num_samples),
    #                   dimnames = list(1:stan_input$time), , 1:num_samples)
    ret$sample <-array(as.numeric(unlist(tmp)),dim=c(nrow(ret$mle[[1]]),ncol(ret$mle[[1]]),num_samples))
    ret$mle <-array(as.numeric(unlist(ret$mle)),dim=c(nrow(ret$mle[[1]]),ncol(ret$mle[[1]]),num_samples))
  }else{
    tmp<-lapply(sammy,function(x,y){
      y - x$sam_x_hat + x$sam_x
    }, y = ret$mle)
    ret$sample <-array(as.numeric(unlist(tmp)),dim=c(nrow(ret$mle), ncol(ret$mle), num_samples))
  }

  class(ret) <- "EnsembleSample"

  rm("KalmanFilter_back")

  return(ret)
}

#'@rdname get_stan_outputs
get_transformed_data <- function(fit){
  stan_input <- fit$stan_input

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
get_parameters <- function(ex.fit, x=1){
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
get_mle <- function(x=1,ex.fit,transformed_data,time) ## get the MLE from the fit
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







#gen_sample <- function(x=1, ex.fit, # variables
#                       bigM, all_eigenvalues_cov, observation_available, time) # constants
#{
#
#  params <- get_parameters(ex.fit, x)
#
#  sam_y <- matrix(0,time,nrow(bigM))
#  sam_x <- matrix(0,time,length(params$x_hat))
#
#  sam_x[1,] <- as.matrix(params$x_hat) + chol(params$SIGMA_init)%*%rnorm(length(params$x_hat))
#  sam_y[1,] <- bigM %*% (sam_x[1,] + params$lt_discrepancies) + rnorm(nrow(bigM),0,sqrt(all_eigenvalues_cov))
#  ch_sigma <- chol(params$SIGMA)
#  for (i in 2:time){
#    sam_x[i,] <- as.matrix(params$AR_params * sam_x[i-1,]) + ch_sigma%*%rnorm(length(params$x_hat))
#    sam_y[i,] <- bigM %*% (sam_x[i,] + params$lt_discrepancies) + rnorm(nrow(bigM),0,sqrt(all_eigenvalues_cov))
#  }
#  sam_x_hat <- KalmanFilter_back(params$AR_params, params$lt_discrepancies,
#                                 all_eigenvalues_cov, params$SIGMA,
#                                 bigM,
#                                 params$SIGMA_init,
#                                 params$x_hat, time, sam_y,
#                                 observation_available)
#  return(list(sam_x=sam_x,sam_x_hat=sam_x_hat))
#}


#get_mle <- function(x=1,ex.fit,bigM,all_eigenvalues_cov,observation_available,time,new_data) ## get the MLE from the fit
#{
#  params <- get_parameters(ex.fit,x)
#  ret <- KalmanFilter_back(params$AR_params, params$lt_discrepancies, all_eigenvalues_cov,
#                           params$SIGMA,bigM,params$SIGMA_init,params$x_hat,time,new_data,observation_available)
#  return(ret)
#}