test_that("test including drivers, with simulated data set",{
  d<- 2
  M <- 2
  MM <- 3
  Time <- 10
  Times_obs <- round(Time * 0.8)
  obs <- matrix(c(0.57,0.75,0.9,1,1.07,1.09,1.04,0.94,0.99,0.54,0.57,1.09,1.67,2.01,2.17,2.33),ncol=d)
  obs.cov <- matrix(c(0.07,-0.02,-0.02,0.01),ncol=d)
  model.cov <- array(0,dim=c(M * MM,d,d))
  models_output <- array(NA,dim=c(M * MM,Time,d))
  models_output[1,,] <- c(-0.14,0.3,0.78,1.19,1.38,1.27,1.04,0.97,1.27,1.85,3.75,3.18,3.01,3.33,3.91,4.41,4.71,4.97,5.42,6.08)
  models_output[2,,] <- c(-0.94,-0.55,-0.25,-0.15,-0.26,-0.49,-0.7,-0.73,-0.52,-0.17,1.48,0.93,0.72,0.98,1.56,2.14,2.54,2.78,2.99,3.25)
  models_output[3,,] <- c(-0.05, 0.35, 0.84, 1.15, 1.18, 1.17, 1.08, 1.16, 1.94, 2.39, 3.93, 3.19, 2.88, 3.17, 3.69, 3.89, 3.99, 4.85, 5.6, 6.2)
  models_output[4,,] <- c(-1.16,-0.75,-0.65,-0.56,-0.66,-0.79,-0.9,-0.43,-0.02, 0.17,1.18,0.27,0.12,0.38,0.46,0.54,0.94,2.08,2.71,3.35)
  models_output[5,,] <- c(-1.14, 1.3, 1.78, 2.19, 2.38, 2.27, 2.04, 1.97, 2.27, 2.85, 4.75, 4.18, 4.01, 4.33, 4.91, 5.41, 5.71, 5.97, 6.42, 7.08)
  models_output[6,,] <- c(0.06, 0.45, 0.75, 0.85, 0.74, 0.51, 0.3, 0.27, 0.48, 0.83, 2.84, 1.93, 1.72, 1.98, 2.56, 3.14, 3.54, 3.78, 3.99, 4.25)

  testplts <- lapply(1:4, function(i) {cbind(reshape2::melt(data.frame(models_output[i,,])), rep(seq(1,10),2), rep(i,20))})
  testdf <- do.call(rbind, testplts)
  names(testdf)[3:4] <- c("year", "combination")
  testdf$combination <- as.factor(testdf$combination)
  testplt <- ggplot2::ggplot(testdf) + geom_line(aes(x = year, y = value, color = variable)) + facet_wrap(vars(combination))
  testplt

  obsdat <- cbind(reshape2::melt(data.frame(obs)), rep(seq(1,8),2))
  names(obsdat)[3] <- "year"
  obsplt <- ggplot2::ggplot(obsdat) + geom_line(aes(x = year, y = value, color = variable))
  obsplt

  model.cov[1,,] <- c(0.06,0,0,0.03)
  model.cov[2,,] <- c(0.04,0.01,0.01,0.02)
  model.cov[3,,] <- c(0.05, 0, 0, 0.02)
  model.cov[4,,] <- c(0.05, -0.03, -0.03, 0.06)
  model.cov[5,,] <- c(0.07,-0.01, -0.01, 0.08)
  model.cov[6,,] <- c(0.06, 0.03, 0.03, 0.05)

  val_obs <- data.frame(obs); cov_obs <- obs.cov
  val_model_11 <- data.frame(models_output[1,,]); cov_model_11 <- model.cov[1,,]
  val_model_21 <- data.frame(models_output[2,,]); cov_model_21 <- model.cov[2,,]
  val_model_12 <- data.frame(models_output[3,,]); cov_model_12 <- model.cov[3,,]
  val_model_22 <- data.frame(models_output[4,,]); cov_model_22 <- model.cov[4,,]
  val_model_13 <- data.frame(models_output[5,,]); cov_model_13 <- model.cov[5,,]
  val_model_23 <- data.frame(models_output[6,,]); cov_model_23 <- model.cov[6,,]

  #Set the dimnames to ensure EcoEnsemble can identify the information.
  SPECIES_NAMES <- c("Species 1", "Species 2")
  dimnames(val_obs) <- list(paste(1:Times_obs), SPECIES_NAMES)
  dimnames(val_model_11) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_21) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_12) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_22) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_13) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_23) <- list(paste(1:Time), SPECIES_NAMES)

  dimnames(cov_obs) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_11) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_21) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_12) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_22) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_13) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_23) <- list(SPECIES_NAMES, SPECIES_NAMES)

  val_model_1 <- list(val_model_11, val_model_12, val_model_13)
  val_model_2 <- list(val_model_21, val_model_22, val_model_23)
  cov_model_1 <- list(cov_model_11, cov_model_12, cov_model_13)
  cov_model_2 <- list(cov_model_21, cov_model_22, cov_model_23)

  for (sampler in c("kalman", "explicit")){
  suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                             simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2", "Driver 3")),
                                                               list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2", "Driver 3"))
                                             ),
                                             priors = EnsemblePrior(d), sampler = sampler,
                                             control = list(adapt_delta = 0.9),chains=1,iter=4,drivers=T))

  priors2 <- EnsemblePrior(2,
                           ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(c(1,2),c(1,1)), cor_params = list(matrix(50, 2, 2),matrix(50, 2, 2))
                                                      , AR_params = c(2, 2)),
                           ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
                           ),
                           sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))
                           ))
  suppressWarnings(fit.opt <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                                 simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2", "Driver 3")),
                                                                   list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2", "Driver 3"))
                                                 ), sampler = sampler,
                                                 priors = priors2,full_sample = F, drivers = T))

  }
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
  MM <- stan_input$MM
  model_num_species <- stan_input$model_num_species
  obs_covariances <- stan_input$obs_covariances
  Ms <- stan_input$Ms
  model_covariances <- stan_input$model_covariances
  observation_times <- stan_input$observation_times
  time <- stan_input$time
  observations <- stan_input$observations
  model_outputs <- stan_input$model_outputs
  n_dri <- stan_input$n_dri
  mod_to_dri <- stan_input$mod_to_dri

  bigM <- matrix(0, N + sum(model_num_species), (M + MM + 2) * N)
  all_eigenvalues_cov <- rep(0, N + sum(model_num_species))
  all_eigenvectors_cov = matrix(0, N + sum(model_num_species), N + sum(model_num_species))
  new_data <- observation_available <- matrix(0, time, N + sum(model_num_species))
  bigM[1:N,1:N] = diag(N)

  tmp_eigen <- eigen(obs_covariances,symmetric=T)
  all_eigenvalues_cov[1:N] <- tmp_eigen$values
  all_eigenvectors_cov[1:N,1:N] <- tmp_eigen$vectors

  bigM[(N + 1):(N + model_num_species[1] ),1:N] <- Ms[1:model_num_species[1],]
  bigM[(N + 1):(N + model_num_species[1] ),(N+1):(2*N)] <-  Ms[1:model_num_species[1],]
  bigM[(N + 1):(N + model_num_species[1] ),(2*N + 1):(3 * N)] <- Ms[1:model_num_species[1],]
  bigM[(N + 1):(N + model_num_species[1]), ((M + 2)*N + (which_fun_subset(1, 1, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(1, 1, mod_to_dri)*N)] <- Ms[1:model_num_species[1],]

  tmp_eigen <- eigen(matrix(model_covariances[1:(model_num_species[1]^2)],model_num_species[1],model_num_species[1]),symmetric=T)
  all_eigenvalues_cov[(N + 1):(N + model_num_species[1] )] <- tmp_eigen$values
  all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] = tmp_eigen$vectors;

  j <- 2

  start_index <- N + sum(model_num_species[1:(j-1)]) + 1
  end_index <-  N + sum(model_num_species[1:j])
  Ms_transformation <- Ms[(start_index - N):(end_index - N),]
  bigM[start_index:end_index, 1:N] <- Ms_transformation
  bigM[start_index:end_index, (N + 1):(2*N)] <- Ms_transformation
  bigM[start_index:end_index, (2*N + 1):(3*N)] <- Ms_transformation
  bigM[start_index:end_index, ((M + 2)*N + (which_fun_subset(1, j, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(1, j, mod_to_dri)*N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:(j-1)]^2) + 1):sum(model_num_species[1:j]^2)],model_num_species[j],model_num_species[j]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

  j <- 3

  start_index <- N + sum(model_num_species[1:(j-1)]) + 1
  end_index <-  N + sum(model_num_species[1:j])
  Ms_transformation <- Ms[(start_index - N):(end_index - N),]
  bigM[start_index:end_index, 1:N] <- Ms_transformation
  bigM[start_index:end_index, (N + 1):(2*N)] <- Ms_transformation
  bigM[start_index:end_index, (2*N + 1):(3*N)] <- Ms_transformation
  bigM[start_index:end_index, ((M + 2)*N + (which_fun_subset(1, j, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(1, j, mod_to_dri)*N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:(j-1)]^2) + 1):sum(model_num_species[1:j]^2)],model_num_species[j],model_num_species[j]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

  i <- 2
  j <- 1

  start_index <-  N + sum(model_num_species[1:sum(n_dri[1:(i-1)])]) + 1
  end_index <- N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + 1)])
  Ms_transformation <- Ms[(start_index - N):(end_index - N),]
  bigM[start_index:end_index,1:N] <- Ms_transformation;
  bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation;
  bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] <- Ms_transformation;
  bigM[start_index:end_index, ((M + 1)*N + (which_fun_subset(i, 1, mod_to_dri)-1)*N + 1):((M + 1)*N + which_fun_subset(i, 1, mod_to_dri)*N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:sum(n_dri[1:(i-1)])]^2) + 1):sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + 1)]^2)], model_num_species[sum(n_dri[1:(i-1)]) + 1], model_num_species[sum(n_dri[1:(i-1)]) + 1]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

  j <- 2

  start_index <-  N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + (j-1))]) + 1
  end_index <- N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + j)])
  Ms_transformation <- Ms[(start_index - N):(end_index - N),]
  bigM[start_index:end_index,1:N] <- Ms_transformation
  bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation
  bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] = Ms_transformation
  bigM[start_index:end_index, ((M + 2)*N + (which_fun_subset(i, j, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(i, j, mod_to_dri)*N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + (j-1))]^2) + 1):sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + j)]^2)], model_num_species[sum(n_dri[1:(i-1)]) + j], model_num_species[sum(n_dri[1:(i-1)]) + j]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

  j <- 3

  start_index <-  N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + (j-1))]) + 1
  end_index <- N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + j)])
  Ms_transformation <- Ms[(start_index - N):(end_index - N),]
  bigM[start_index:end_index,1:N] <- Ms_transformation
  bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation
  bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] = Ms_transformation
  bigM[start_index:end_index, ((M + 2)*N + (which_fun_subset(i, j, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(i, j, mod_to_dri)*N)] <- Ms_transformation

  tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + (j-1))]^2) + 1):sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + j)]^2)], model_num_species[sum(n_dri[1:(i-1)]) + j], model_num_species[sum(n_dri[1:(i-1)]) + j]),symmetric=T)
  all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
  all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors

  for (i in 1:time){
    observation_available[i,1:N] <- rep(observation_times[i,1],N)
    if (M > 0){
      observation_available[i,(N+1):(N + model_num_species[1])] <- rep(observation_times[i,2],model_num_species[1])
      for (j in 2:M){
        observation_available[i,(N + sum(model_num_species[1:sum(n_dri[1:(j-1)])]) + 1):(N + sum(model_num_species[1:sum(n_dri[1:j])]))] <- rep(observation_times[i,j + 1],sum(model_num_species[(sum(n_dri[1:(j-1)]) + 1):sum(n_dri[1:j])]))
      }
    }
    new_data[i,] <- t(all_eigenvectors_cov) %*% matrix(c(observations[i,],model_outputs[i,]),ncol=1)
  }

  bigM <- t(all_eigenvectors_cov) %*% bigM
  ret <- list(bigM=bigM,
              all_eigenvalues_cov=all_eigenvalues_cov,
              new_data=new_data,
              observation_available=observation_available)
  expect_equal(get_transformed_data_dri(fit),ret)

  #test get_parameters
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
  transformed_data <- get_transformed_data_dri(fit)
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
  transformed_data <- get_transformed_data_dri(fit.opt)
  params <- get_parameters(fit.opt@point_estimate$par)

  ret1 <- KalmanFilter_back(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                            params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                            transformed_data$new_data, transformed_data$observation_available)

  expect_equal(back_kalmanFilter_R(params$AR_params, params$lt_discrepancies, transformed_data$all_eigenvalues_cov,
                                   params$SIGMA, transformed_data$bigM, params$SIGMA_init, params$x_hat,time,
                                   transformed_data$new_data, transformed_data$observation_available),ret1)
})
