validate_data_dri <- function(observations, simulators, priors, MMod){
  validate_observations(observations)

  if(!is.list(simulators)){
    stop(paste0("Simulator data should be passed through as a list, with each element also a list corresponding to each model. For each sublist corresponding to a model, the first element should be a further list containing the data for each of the allowed driver combinations with the model, and the second element a list containing the corresponding covariances. The third element element is the optional model name and the fourth element is a list of optional driver names (for each allowed combination)."))
  }

  if (!missing(MMod)){
    MM <- max(unique(unlist(MMod)))
    n_dri <- sapply(MMod, length)
  }
  else {
    MM <- length(simulators[[1]][[1]])
    MMod <- replicate(M, 1:MM, simplify = F)
    n_dri <- rep(MM, M)
  }

  for(i in 1:M){
    simulator <- simulators[[i]]
    sim_name <- paste0("One of your simulators [[", i, "]]")
    if(length(simulator) >= 3){
      sim_name <- simulator[[3]]
    }
    if(!is.list(simulator)){
      stop(paste0(sim_name, "is not a list. Each simulator should be passed as a list, the first element of which should be a further list containing the data for each of the allowed driver combinations with the simulator, and the second element should be a list containing the corresponding covariances. The (optional) third element element is the simulator name and the (optional) fourth element is a list driver names (for each allowed combination with the given simulator)."))
    }
    #Check supplied model/driver combinations are compatible
    if (!missing(MMod)){
      if (length(simulator[[1]]) != n_dri[i] || length(simulator[[2]]) != n_dri[[i]]){
        stop(paste0(sim_name), "has a number of driver combinations incompatible with user-supplied number of combinations:", n_dri[i])
      }
    }
    else {
      if (length(simulator[[1]]) != MM || length(simulator[[2]]) != MM){
        stop(paste0(sim_name), "does not have data for all driven combinations. Data for all", MM,"drivers should be included.")
      }
    }
    for (j in 1:n_dri[i]){
      driver <- simulator[[1]][[j]]
      driver_sigma <- simulator[[2]][[j]]
      sim_name <- paste0("One of your simulator/driver combinations [[", i, "]][[", j,"]]")
      dri_name <- paste0("[[", j,"]]")
      if(length(simulator) >= 3){
        sim_name <- simulator[[3]]
        dri_name <- paste0(" + ","driver [[", j,"]]")
        if (length(simulator) = 4){
          dri_name <- paste0(" + ", simulator[[4]][[j]])
        }
        sim_name <- paste0(sim_name, dri_name)
      }
      sim_name <- paste0(sim_name, dri_name)
      validate_simulator(list(driver, driver_sigma), colnames(observations[[1]]), sim_name)
    }
  }
}
