generate_simulator_stan_data_dri <- function(observations, simulators, MMod){
  N <- ncol(observations[[1]])
  M <- length(simulators)
  if (!missing(MMod)){
    MM <- max(unique(unlist(MMod)))
    n_dri <- sapply(MMod, length)
  }
  else {
    MM <- length(simulators[[1]][[1]])
    MMod <- replicate(M, 1:MM, simplify = F)
    n_dri <- rep(MM, M)
  }
  times <- rownames(observations[[1]])
  model_num_species <- c()
  Ms <- matrix(NA, nrow = 0, ncol = N)
  model_covariances <- c()
  observed_species <- colnames(observations[[1]])
  mod_to_dri <- matrix(0, nrow = M, ncol = MM)
  for (i in 1:M) {
    sim <- simulators[[i]]
    for (j in 1:n_dri[i]) {
      model_output <- sim[[1]][[j]]
      times <- unique(append(times, rownames(model_output)))
      mod_num_spec <- ncol(model_output)
      model_num_species <- c(model_num_species, mod_num_spec)
      sim_species <- colnames(model_output)
      Mij <- matrix(0, nrow = mod_num_spec, ncol = N)
      for (k in 1:mod_num_spec){
        for (l in 1:N){
          Mij[k,l] <- sim_species[k] == observed_species[l]
        }
      }
      Ms <- rbind(Ms, Mij)
      ## covariances
      if (is.matrix(sim[[2]][[j]])){
        model_covariances <- c(model_covariances, as.vector(sim[[2]][[j]]))
      }else{
        model_covariances <- c(model_covariances, as.numeric(sim[[2]][[j]]))
      }
      dri_index <- MMod[[i]][j]
      mod_to_dri[i, dri_index] <- dri_index
    }
  }

  times <- sort(as.integer(times))

  #Things that can't be done on the first pass
  observation_times <- matrix(NA, nrow=length(times), ncol=0)
  observation_times <- cbind(observation_times, times %in% rownames(observations[[1]]))
  for (i in 1:M) {
    sim <- simulators[[i]]
    model_output <- sim[[1]][[1]]
    present_data <- as.numeric(times %in% rownames(model_output))
    observation_times <- cbind(observation_times, present_data)
  }

  model_outputs <- matrix(NA, nrow=length(times), ncol=0)
  for (i in 1:M) {
    sim <- simulators[[i]]
    for (j in 1:n_dri[i]){
      model_output <- sim[[1]][[j]]
      y_i <- matrix(0, nrow = length(times), ncol=ncol(model_output))
      for (k in 1:length(times)){
        year = times[k]
        dat_for_year <- model_output[year == rownames(model_output), ]
        if (nrow(dat_for_year) == 1 || length(dat_for_year) == 1){
          y_i[k, ] <- unlist(dat_for_year)
        }
      }
      model_outputs <- cbind(model_outputs, y_i)
    }
  }

  obs_data <- observations[[1]]
  obs_data_all <- matrix(0, nrow = length(times), ncol=N)
  for (k in 1:length(times)){
    year = times[k]
    obs_for_year <- obs_data[year == rownames(obs_data), ]
    if (nrow(obs_for_year) == 1 || length(obs_for_year) == 1){
      obs_data_all[k, ] <- unlist(obs_for_year)
    }
  }
  obs_covariances <- observations[[2]]

  return(list(N = N,
              time = length(times),
              M = M,
              MM=MM,
              n_dri = n_dri,
              mod_to_dri = mod_to_dri,
              model_num_species = model_num_species,
              Ms = Ms,
              observation_times = observation_times,
              model_outputs = model_outputs,
              model_covariances = model_covariances,
              observations = obs_data_all,
              obs_covariances = obs_covariances))
}
