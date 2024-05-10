# These are the currently accepted priors for correlation matrices. If you add anything here, you will also need to change the validate_parametrisation_form function in validation_prior_functions.R
CORRELATIONS_PRIOR_LKJ <- 0
CORRELATIONS_PRIOR_INV_WISHART <- 1
CORRELATIONS_PRIOR_BETA <- 2
CORRELATIONS_PRIOR_HIERARCHICAL <- 3
CORRELATIONS_PRIOR_BETA_CONJUGATE <- 4

correlation_form_prior <- function(prior_choice){
  return(
    switch(EXPR = tolower(prior_choice),
           "lkj" = CORRELATIONS_PRIOR_LKJ,
           "inv_wishart" = CORRELATIONS_PRIOR_INV_WISHART,
           "beta" = CORRELATIONS_PRIOR_BETA,
           "hierarchical" = CORRELATIONS_PRIOR_HIERARCHICAL,
           "hierarchical_beta_conjugate" = CORRELATIONS_PRIOR_BETA_CONJUGATE)
  )
}


correlation_prior <- function(params, form, type, allow_hierarchical){
  retv <- rep(NA, 5)
  names(retv) <- c(paste0(type, "_lkj"), paste0(type, "_wish_nu"),
                   paste0(type, "_wish_sigma"),paste0(type, "_beta_1"),
                   paste0(type, "_beta_2"))


  if(allow_hierarchical){
    retv <- append(retv, NA)
    names(retv)[6] <- paste0(type, "_hierarchical_beta_hyper_params")
  }

  ret <- as.list(retv)

  ret[[paste0(type, "_lkj")]] <- numeric(0)
  ret[[paste0(type, "_wish_nu")]] <- numeric(0)
  ret[[paste0(type, "_wish_sigma")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_1")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_2")]] <- matrix(0,0,0)

  if(allow_hierarchical){
    ret[[paste0(type, "_hierarchical_beta_hyper_params")]] <- numeric(0)
  }


  if (form == CORRELATIONS_PRIOR_LKJ){
    ret[[paste0(type, "_lkj")]] <- array(params, 1)
  }else if(form == CORRELATIONS_PRIOR_INV_WISHART){
    ret[[paste0(type, "_wish_nu")]] <-array(params[[1]], 1)
    ret[[paste0(type, "_wish_sigma")]] <- params[[2]]
  }else if(form == CORRELATIONS_PRIOR_BETA){
    ret[[paste0(type, "_beta_1")]] <- params[[1]]
    ret[[paste0(type, "_beta_2")]] <- params[[2]]
  }else if(form == CORRELATIONS_PRIOR_HIERARCHICAL
           && allow_hierarchical){
    ret[[paste0(type, "_hierarchical_beta_hyper_params")]] <- unlist(params)
  }else if(form == CORRELATIONS_PRIOR_BETA_CONJUGATE
           && allow_hierarchical){
    ret[[paste0(type, "_hierarchical_beta_hyper_params")]] <- unlist(params)
  }
  return(ret)
}

#If the user only specifies one variance parameter common for each species, then we should repeat this for each species so that Stan can understand.
repeat_variance_priors <- function(d, priors){

  prior_var_params <- priors@var_params

  ret <- prior_var_params

  if(length(prior_var_params[[1]]) == 1 && length(prior_var_params[[2]] == 1)){
    if(d == 1){
      ret[[1]] <- array(prior_var_params[[1]], dim = 1)
      ret[[2]] <- array(prior_var_params[[2]], dim = 1)
    }else{
      ret[[1]] <- rep(prior_var_params[[1]], d)
      ret[[2]] <- rep(prior_var_params[[2]], d)
    }
  }

  return(ret)
}

