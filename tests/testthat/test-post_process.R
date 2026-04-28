test_that("test with simulated data set",{
  d<- 2
  M <- 2
  Time <- 10
  Times_obs <- round(Time * 0.8)
  obs <- matrix(c(0.57,0.75,0.9,1,1.07,1.09,1.04,0.94,0.99,0.54,0.57,1.09,1.67,2.01,2.17,2.33),ncol=d)
  obs.cov <- matrix(c(0.07,-0.02,-0.02,0.01),ncol=d)
  model.cov <- array(0,dim=c(M,d,d))
  models_output <- array(NA,dim=c(M,Time,d))
  models_output[1,,] <- c(-0.14,0.3,0.78,1.19,1.38,1.27,1.04,0.97,1.27,1.85,3.75,3.18,3.01,3.33,3.91,4.41,4.71,4.97,5.42,6.08)
  models_output[2,,] <- c(-0.94,-0.55,-0.25,-0.15,-0.26,-0.49,-0.7,-0.73,-0.52,-0.17,1.48,0.93,0.72,0.98,1.56,2.14,2.54,2.78,2.99,3.25)
  model.cov[1,,] <- c(0.06,0,0,0.03)
  model.cov[2,,] <- c(0.04,0.01,0.01,0.02)

  val_obs <- data.frame(obs); cov_obs <- obs.cov
  val_model_1 <- data.frame(models_output[1,,]); cov_model_1 <- model.cov[1,,]
  val_model_2 <- data.frame(models_output[2,,]); cov_model_2 <- model.cov[2,,]


  #Set the dimnames to ensure EcoEnsemble can identify the information.
  SPECIES_NAMES <- c("Species 1", "Species 2")
  dimnames(val_obs) <- list(paste(1:Times_obs), SPECIES_NAMES)
  dimnames(val_model_1) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_2) <- list(paste(1:Time), SPECIES_NAMES)


  dimnames(cov_obs) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_1) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_2) <- list(SPECIES_NAMES, SPECIES_NAMES)

  for (sampler in c("kalman", "explicit")){
  suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                           simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                             list(val_model_2, cov_model_2, "Model 2")
                                           ), sampler = sampler,
                                           priors = EnsemblePrior(d),
                                           control = list(adapt_delta = 0.9),chains=1,iter=4))
  priors2 <- EnsemblePrior(2,
	ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(c(1,2),c(1,1)), cor_params = list(matrix(50, 2, 2),matrix(50, 2, 2))
  , AR_params = c(2, 2)),
	ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
  ),
	sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))
  ))
  suppressWarnings(fit.opt <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                                 list(val_model_2, cov_model_2, "Model 2")
                                               ), sampler = sampler,
                                               priors = priors2,full_sample = F))


  samples <- generate_sample(fit)
  expect_true(ggplot2::is_ggplot(plot(samples,quantiles = c(0.05, 0.95))))
  samples.opt <- generate_sample(fit.opt,num_samples = 500)
  expect_true(ggplot2::is_ggplot(plot(samples.opt,quantiles = c(0.05, 0.95))))
  err_mess <- paste0("Invalid variable. This should be the name of a variable or an index less than ",
                     d + 1)
  expect_error(plot(samples.opt,quantiles = c(0.05, 0.95),variable = 4),err_mess)
  ### test get_transformed_data

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

  bigM[(N + 1):(N + model_num_species[1] ),1:N] <- Ms[1:model_num_species[1],]
  bigM[(N + 1):(N + model_num_species[1] ),(N+1):(2*N)] <-  Ms[1:model_num_species[1],]
  bigM[(N + 1):(N + model_num_species[1] ),(2*N + 1):(3 * N)] <- Ms[1:model_num_species[1],]

  tmp_eigen <- eigen(matrix(model_covariances[1:(model_num_species[1]^2)],model_num_species[1],model_num_species[1]),symmetric=T)
  all_eigenvalues_cov[(N + 1):(N + model_num_species[1] )] <- tmp_eigen$values
  all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] = tmp_eigen$vectors;

  i <- 2

  start_index <-  N + sum(model_num_species[1:(i-1)]) + 1
  end_index <- N + sum(model_num_species[1:i])
  Ms_transformation <- Ms[(sum(model_num_species[1:(i-1)]) + 1):sum(model_num_species[1:i]),]
  bigM[start_index:end_index,1:N] <- Ms_transformation
  bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation
  bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(cumsum(model_num_species^2)[i-1] + 1):(cumsum(model_num_species^2)[i])],model_num_species[i],model_num_species[i]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

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
  ret <- list(bigM=bigM,
              all_eigenvalues_cov=all_eigenvalues_cov,
              new_data=new_data,
              observation_available=observation_available)
  expect_equal(get_transformed_data(fit),ret)


  ex.fit <- rstan::extract(fit@samples)
  x <- 1
  ret <- list()
  ret$x_hat <- ex.fit$x_hat[x,]
  ret$SIGMA_init <- ex.fit$SIGMA_init[x,,]
  ret$AR_params <- ex.fit$AR_params[x,]
  ret$new_x <- ex.fit$x_hat[x,]
  ret$SIGMA <- ex.fit$SIGMA[x,,]
  ret$lt_discrepancies <- ex.fit$lt_discrepancies[x,]
  expect_equal(get_parameters(ex.fit),ret)

  ### test get_mle
  ex.fit <- rstan::extract(fit@samples)
  transformed_data <- get_transformed_data(fit)
  time <- fit@ensemble_data@stan_input$time

  params <- get_parameters(ex.fit,x=1)

  ret <- KalmanFilter_back(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                           params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                           transformed_data$new_data, transformed_data$observation_available)

  expect_equal(get_mle(x=1, ex.fit, transformed_data = transformed_data,
                       time = time),ret)
  ### test kalman filer
  back_kalmanFilter_R <- function(rhos, dee, R, Q, C, P, xhat, Time, y, obs){
    ## initials
    A <- as.matrix(rhos) %*% t(as.matrix(rhos))
    Gs <- array(0,dim=c(length(rhos),ncol(y),Time))
    Ps <- array(0,dim=c(length(rhos),length(rhos),Time))
    xhats <- matrix(0,length(rhos),Time)
    er <- Qinv <- matrix(0,Time,ncol(y))

    for (i in 1:Time)
    {
      P <-  P * A + Q
      xhat <- rhos * xhat
      Ps[,,i] <- P
      xhats[,i] <- xhat
      I <- diag(length(xhat))
      for (j in 1:ncol(y))
      {
        if(obs[i,j]==1.0){
          g <- P %*%  C[j,]
          Qinv[i,j] <- as.numeric(1/(t(C[j,]) %*% g + R[j]))
          er[i,j] <- y[i,j] - sum((xhat + dee) * C[j,])
          Gs[,j,i] <- Qinv[i,j] * g
          xhat <- xhat + Gs[,j,i] * er[i,j]
          P <- P - Gs[,j,i] %*% t(g)
        }
      }
    }
    ret <- matrix(0,Time,length(xhat))
    r <- rep(0,length(rhos))
    for (i in Time:1){
      for(j in ncol(y):1){
        if(obs[i,j]==1.0){
          L <- I - Gs[,j,i] %*% t(C[j,])
          r <- C[j,] * Qinv[i,j] * er[i,j] + t(L) %*% r
        }
      }
      ret[i,] <- xhats[,i] + t(Ps[,,i] %*% r)
      r <- rhos * as.numeric(r)
    }
    return(ret)
  }


  ### test kalman filter test
  time <- fit.opt@ensemble_data@stan_input$time
  transformed_data <- get_transformed_data(fit.opt)
  params <- get_parameters(fit.opt@point_estimate$par)

  ret1 <- KalmanFilter_back(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                            params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                            transformed_data$new_data, transformed_data$observation_available)

  expect_equal(back_kalmanFilter_R(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                                   params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                                   transformed_data$new_data, transformed_data$observation_available),ret1)
}})
