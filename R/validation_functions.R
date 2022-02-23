validate_compatibility <- function(name, pair){
  covar <- as.matrix(pair[[2]])

  if (!is.square.matrix(covar)
      || !is.positive.definite(covar)){
    msg <- paste0(name, " covariance matrix is not a positive definite square matrix.")
    stop(msg)
  }

  df <- as.matrix(pair[[1]])
  if (ncol(df) != ncol(covar)){
    msg <- paste0(name, " data has ", ncol(df), "columns but the covariance matrix has ", ncol(covar), "columns. Variables for the covariance matrix should match the observed data.")
    stop(msg)
  }


  data_variables <- colnames(df)
  cova_variables <- colnames(covar)

  if (length(data_variables) == 0){

    msg <- paste0(name, " data frame variables are not named. Each column of the data frame and the covariance matrix should be named appropriately and match each other.")
    stop(msg)
  }

  if (length(cova_variables) == 0){
    msg <- paste(name ," covariance matrix variables are not",
                 " named. Each column of the data frame",
                 " and the covariance matrix should be named",
                 " appropriately and match each other.", sep="")
    stop(msg)
  }

  if (!all(cova_variables == data_variables)){
    #TODO Split nicely over multiple lines???
    msg <- paste(paste(name, " covariance matrix variables do not",
                       " match the associated data frame. Each column",
                       " of the data frame",
                       " and covariance matrix should match each other",
                       " and be in the same order.", sep=""),
                 "Data variables:",
                 paste(data_variables, sep="", collapse=', '),
                 "Covariance matrix variables: ",
                 paste(cova_variables, sep="", collapse=', '), sep =" ")
    stop(msg)
  }
}

validate_input_data <- function(name, input_data){
  if(!is.list(input_data)||
     (length(input_data) != 2 && length(input_data) != 3) ||
     !(is.matrix(input_data[[1]]) || is.data.frame(input_data[[1]])) ||
     !(is.matrix(input_data[[2]]) || is.data.frame(input_data[[2]]))){
    msg <- paste0(name, " data is not in the correct form. This data should be ",
                  "passed through as a list, the first element being a matrix or data ",
                  "frame containing the data, the second element being the covariance ",
                  "matrix as a matrix or data frame. The name can also be passed through",
                  "as an optional 3rd element.")
    stop(msg)

  }
}

validate_observations <- function(observations){
  validate_input_data("Observation", observations)
  validate_compatibility("Observation", observations)
}

validate_simulator <- function(simulator, observation_variables, name){
  #Check that the data is at least a list
  validate_input_data(name, simulator)
  # Check that it's internally compatible
  validate_compatibility(name, simulator)

  # Check also that it's compatible with the observations
  df <- simulator[[1]]
  if (!all(colnames(df) %in% observation_variables)){
    #TODO: Wrap over a few lines.
    msg <- paste(name, "contains variables that are not present in the observed data. Model outputs",
                 " should match the observed data. Model variables:",
                 paste(colnames(df), sep="", collapse=', '), ". Observation variables:",
                 paste(observation_variables, sep="", collapse=', '), sep="")
    stop(msg)
  }

}


validate_correlation_priors <- function(type_str, form, params, N) {
  form <- tolower(form)
  if (!(form %in% c("lkj", "inv_wishart", "beta"))){
    msg <- paste0("Invalid prior choice on ", type_str," correlation matrix.",
                  " Correlation priors should be one of 'lkj', 'inv_wishart',",
                  " or 'beta'. Prior choice: ", form)
    stop(msg)
  }

  msg <- ""
  if(form == "lkj" &
     length(params) != 1){
    msg <- paste0("Invalid lkj parameter for ", type_str, " correlation matrix prior.",
                  " This should be a list with exactly one real number.")

  }else if (form == "inv_wishart" &
            (length(params) != 2 ||
              !is.matrix(params[[2]] ))){
    #TODO: Come back to this message.
    msg <- paste0("Invalid inverse Wishart parameters for ", type_str, " correlation",
                  " matrix prior.")
  }else if (form == "beta" &
            (length(params) !=2 ||
             !is.matrix(params[[1]]) ||
             !is.matrix(params[[2]]) ||
             !is.square.matrix(params[[1]]) ||
             !is.square.matrix(params[[2]]) ||
             !is.symmetric.matrix(params[[1]]) || #TODO:What are the exact conditions we should check here? Do we need symmetry?
             !is.symmetric.matrix(params[[2]]) ||
             !(dim(params[[1]])[1] == N) ||
             !(dim(params[[2]])[1] == N))){
    msg <- paste0("Invalid beta parameters for ", type_str, " correlation matrix priors.",
                  " These should be square, symmetric matrices with the same dimension as the",
                  " number of model outputs.")
  }
  if(msg != ""){
    stop(msg)
  }

}

validate_priors <- function(N, priors){

  #Correlation priors
  cor_pri <- priors$prior_correlations
  validate_correlation_priors("individual short-term", cor_pri$ind_st_cor_form, cor_pri$ind_st_cor_params, N)
  validate_correlation_priors("individual long-term",  cor_pri$ind_lt_cor_form, cor_pri$ind_lt_cor_params, N)
  validate_correlation_priors("shared short-term",     cor_pri$sha_st_cor_form, cor_pri$sha_st_cor_params, N)

  #Other priors
  if(length(priors$prior_ind_lt_var_a) != N ||
     length(priors$prior_ind_lt_var_b) != N){
    msg <- paste0("Invalid priors for the variance of the individual long-term",
                  " discrepancies. These should be vectors of length ", N)
    stop(msg)

  }
  #TODO: Should we be doing validation on the following?
  #priors$prior_ind_st_var_a
  #priors$prior_ind_st_var_b
  #priors$prior_sha_st_var_exp
}

validate_data <- function(observations, simulators, priors){
  validate_observations(observations)

  if(!is.list(simulators)){
   stop(paste0("Simulator data should be passed through as a list, with each ",
               "element containing a list of the data for a given simulator."))
  }

  #for(simulator in simulators){
  for(i in 1:length(simulators)){
    simulator <- simulators[[i]]
    sim_name <- paste0("One of your simulator's [[", i, "]]")
    if(length(simulator) == 3){
      sim_name <- simulator[[3]]
    }
    validate_simulator(simulator, colnames(observations[[1]]), sim_name)
  }
  validate_priors(ncol(observations[[1]]), priors)
}
