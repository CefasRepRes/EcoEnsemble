which_fun_subset <- function(i, j, mod_to_esm){
  ret <- 0
  num <- 0
  while (num < j){
    ret <- ret + 1
    if (mod_to_esm[i,ret] == ret){
      num <- num + 1
    }
  }
  return(ret)
}


#'@rdname get_stan_outputs
#'@export
get_transformed_data_dri <- function(fit){
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
  all_eigenvalues_cov[1:N] = tmp_eigen$values
  all_eigenvectors_cov[1:N,1:N] = tmp_eigen$vectors
  if (M > 0) {
    if (MM > 0) {
      bigM[(N + 1):(N + model_num_species[1]), 1:N] <- Ms[1:model_num_species[1],]
      bigM[(N + 1):(N + model_num_species[1]), (N + 1):(2*N)] <- Ms[1:model_num_species[1],]
      bigM[(N + 1):(N + model_num_species[1]), (2*N + 1):(3*N)] <- Ms[1:model_num_species[1],]
      bigM[(N + 1):(N + model_num_species[1]), ((M + 2)*N + (which_fun_subset(1, 1, mod_to_dri)-1)*N + 1):((M + 2)*N + which_fun_subset(1, 1, mod_to_dri)*N)] <- Ms[1:model_num_species[1],]

      tmp_eigen <- eigen(matrix(model_covariances[1:(model_num_species[1]^2)],model_num_species[1],model_num_species[1]),symmetric=T)
      all_eigenvalues_cov[(N + 1):(N + model_num_species[1])] <- tmp_eigen$values
      all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] <- tmp_eigen$vectors

      if (MM > 1) {
        for (j in 2:n_dri[1]) {
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
        }
      }
    }
    else {
      bigM[(N + 1):(N + model_num_species[1]), 1:N] <- Ms[1:model_num_species[1],]
      bigM[(N + 1):(N + model_num_species[1]), (N + 1):(2*N)] <- Ms[1:model_num_species[1],]
      bigM[(N + 1):(N + model_num_species[1]), (2*N + 1):(3*N)] <- Ms[1:model_num_species[1],]
      tmp_eigen <- eigen(matrix(model_covariances[1:(model_num_species[1]^2)],model_num_species[1],model_num_species[1]),symmetric=T)
      all_eigenvalues_cov[(N + 1):(N + model_num_species[1])] <- tmp_eigen$values
      all_eigenvectors_cov[(N + 1):(N + model_num_species[1]),(N + 1):(N + model_num_species[1])] <- tmp_eigen$vectors
    }
    if (M > 1) {
      if (MM > 0) {
        for (i in 2:M) {
          start_index <-  N + sum(model_num_species[1:sum(n_dri[1:(i-1)])]) + 1
          end_index <- N + sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + 1)])
          Ms_transformation <- Ms[(start_index - N):(end_index - N),]
          bigM[start_index:end_index,1:N] <- Ms_transformation;
          bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation;
          bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] <- Ms_transformation
          bigM[start_index:end_index, ((M + 1)*N + (which_fun_subset(i, 1, mod_to_dri)-1)*N + 1):((M + 1)*N + which_fun_subset(i, 1, mod_to_dri)*N)] <- Ms_transformation
          tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:sum(n_dri[1:(i-1)])]^2) + 1):sum(model_num_species[1:(sum(n_dri[1:(i-1)]) + 1)]^2)], model_num_species[sum(n_dri[1:(i-1)]) + 1], model_num_species[sum(n_dri[1:(i-1)]) + 1]),symmetric=T)
          all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
          all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors
        }
        if (MM > 1) {
          for (i in 2:M) {
            for (j in 2:n_dri[i]) {
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
            }
          }
        }
      }
      else {
        for (i in 2:M) {
          start_index <-  N + sum(model_num_species[1:(i-1)]) + 1
          end_index <- N + sum(model_num_species[1:i])
          Ms_transformation <- Ms[(sum(model_num_species[1:(i-1)]) + 1):sum(model_num_species[1:i]),]
          bigM[start_index:end_index,1:N] <- Ms_transformation
          bigM[start_index:end_index,(N+1):(2*N)] <- Ms_transformation
          bigM[start_index:end_index,((i + 1)*N + 1):((i + 2) * N)] <- Ms_transformation
          tmp_eigen <- eigen(matrix(model_covariances[(sum(model_num_species[1:(i-1)]^2) + 1):sum(model_num_species[1:i]^2)],model_num_species[i],model_num_species[i]),symmetric=T)
          all_eigenvalues_cov[start_index:end_index] <- tmp_eigen$values
          all_eigenvectors_cov[start_index:end_index, start_index:end_index] <- tmp_eigen$vectors
        }
      }
    }
  }

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


  bigM = t(all_eigenvectors_cov) %*% bigM;

  return(
    list(bigM=bigM,
         all_eigenvalues_cov=all_eigenvalues_cov,
         new_data=new_data,
         observation_available=observation_available)
  )
}

