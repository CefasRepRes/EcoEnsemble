generate_simulator_stan_data <- function(observations, simulators){
  N <- ncol(observations[[1]])
  M <- length(simulators)
  times <- rownames(observations[[1]])
  model_num_species <- rep(NA, M)
  Ms  <- matrix(NA, nrow=0, ncol=N)
  model_covariances <- c()
  for (i in 1:M) {
    sim <- simulators[[i]]
    model_ouput <- sim[[1]]
    times <- unique(append(times, rownames(model_ouput)))
    model_num_species[i] <- ncol(model_ouput)

    Mi <- matrix(0, nrow=model_num_species[i], ncol=N)
    observed_species <- colnames(observations[[1]])
    sim_species      <- colnames(model_ouput)
    for(k in 1:model_num_species[i]){
      for(l in 1:N){
        Mi[k,l] <- sim_species[k] == observed_species[l]
      }
    }
    Ms <- rbind(Ms, Mi)

    model_covariances <- append(model_covariances, as.numeric(sim[[2]]))

  }
  times <- sort(as.integer(times))

  #Things that can't be done on the first pass
  observation_times <- matrix(NA, nrow=length(times), ncol=0)
  model_outputs <- matrix(NA, nrow=length(times), ncol=0)
  observation_times <- cbind(observation_times, times %in% rownames(observations[[1]]))
  for (i in 1:M) {
    sim <- simulators[[i]]
    model_ouput <- sim[[1]]
    present_data <- as.numeric(times %in% rownames(model_ouput))
    observation_times <- cbind(observation_times, present_data )


    y_i<- matrix(0, nrow = length(times), ncol=ncol(model_ouput))
    for(k in 1:length(times)){
      year = times[k]
      dat_for_year <- model_ouput[year == rownames(model_ouput), ]
      if(nrow(dat_for_year) == 1 ||
         length(dat_for_year) == 1){
        y_i[k, ] <- unlist(dat_for_year)
      }
    }
    model_outputs <- cbind(model_outputs, y_i)

  }

  if(M == 1){
    model_num_species = array(model_num_species, 1)
  }

  obs_data <- observations[[1]]
  obs_data_all <- matrix(0, nrow = length(times), ncol=N)
  for(k in 1:length(times)){
    year = times[k]
    obs_for_year <- obs_data[year == rownames(obs_data), ]
    if(nrow(obs_for_year) == 1 ||
       length(obs_for_year) == 1){
      obs_data_all[k, ] <- unlist(obs_for_year)
    }
  }
  obs_covariances <- observations[[2]]
  return(list(N = N,
              time = length(times),
              M = M,
              model_num_species = model_num_species,
              Ms = Ms,
              observation_times = observation_times,
              model_outputs = model_outputs,
              model_covariances = model_covariances,
              observations = obs_data_all,
              obs_covariances = obs_covariances))
  }




