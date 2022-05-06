# These are the currently accepted priors for correlation matrices. If you add
# anything here, you will also need to change the validate_correlation_priors
# function in validation_functions.R
CORRELATIONS_PRIOR_LKJ <- 0
CORRELATIONS_PRIOR_INV_WISHART <- 1
CORRELATIONS_PRIOR_BETA <- 2
CORRELATIONS_PRIOR_HIERARCHICAL <- 3

correlation_form_prior <- function(prior_choice){
  return(
    switch(EXPR = tolower(prior_choice),
           "lkj" = CORRELATIONS_PRIOR_LKJ,
           "inv_wishart" = CORRELATIONS_PRIOR_INV_WISHART,
           "beta" = CORRELATIONS_PRIOR_BETA,
           "hierarchical" = CORRELATIONS_PRIOR_HIERARCHICAL)
  )
}

correlation_prior <- function(params, form, type){
  retv <- rep(NA, 5)
  names(retv) <- c(paste0(type, "_lkj"), paste0(type, "_wish_nu"),
                   paste0(type, "_wish_sigma"),paste0(type, "_beta_1"),
                   paste0(type, "_beta_2"))

  #JM 06/05/2022: Now allow for hierarchical priors, but only for the individual short-terms
  #TODO: Shouldn't really have magic sneaky things for different types...
  is_individual_short_term = (type == "prior_ind_st_cor")
  if(is_individual_short_term){
    retv <- rep(NA, 6)
    names(retv) <- c(paste0(type, "_lkj"), paste0(type, "_wish_nu"),
                     paste0(type, "_wish_sigma"),paste0(type, "_beta_1"),
                     paste0(type, "_beta_2"),
                     paste0(type, "_hierarchical_beta_hyper_params"))
  }

  ret <- as.list(retv)

  ret[[paste0(type, "_lkj")]] <- numeric(0)
  ret[[paste0(type, "_wish_nu")]] <- numeric(0)
  ret[[paste0(type, "_wish_sigma")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_1")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_2")]] <- matrix(0,0,0)

  if(is_individual_short_term){
    ret[[paste0(type, "_hierarchical_beta_hyper_params")]] <- numeric(4)
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
           && is_individual_short_term){
    #JM 06/05/2022
    ret[[paste0(type, "_hierarchical_beta_hyper_params")]] <- params
  }else{
    stop("Invalid correlation parameters! Did you try hierarchical parameters for something other than the individual short-term?")
  }
  return(ret)
}


generate_correlation_priors_stan_data <- function(prior_correlations){
  # Create the prior form identifier e.g. "lkj" becomes 0
  form_prior_ind_st <- correlation_form_prior(prior_correlations$ind_st_cor_form)
  form_prior_ind_lt <- correlation_form_prior(prior_correlations$ind_lt_cor_form)
  form_prior_sha_st <- correlation_form_prior(prior_correlations$sha_st_cor_form)

  priors_data <- list(
    form_prior_ind_st = form_prior_ind_st,
    form_prior_ind_lt = form_prior_ind_lt,
    form_prior_sha_st = form_prior_sha_st
  )

  priors_data <- append(priors_data, correlation_prior(prior_correlations$ind_st_cor_params, form_prior_ind_st, "prior_ind_st_cor"))
  priors_data <- append(priors_data, correlation_prior(prior_correlations$ind_lt_cor_params, form_prior_ind_lt, "prior_ind_lt_cor"))
  priors_data <- append(priors_data, correlation_prior(prior_correlations$sha_st_cor_params, form_prior_sha_st, "prior_sha_st_cor"))

  return(priors_data)
}



# If only one of each inverse-gamma variable is provided, then repeat it for every species considered.
repeat_priors <- function(d, prior_params, nm, nm_var){
    if (is.list(prior_params) && length(prior_params) == 3 && is.list(prior_params[[2]]) && length(prior_params[[2]]) == 2){
      if(length(prior_params[[2]][[1]]) == 1 && length(prior_params[[2]][[2]]) == 1){
        prior_params[[2]][[1]] <- rep(prior_params[[2]][[1]], d)
        prior_params[[2]][[2]] <- rep(prior_params[[2]][[2]], d)
      }

      #JM: 06-05-2022 Short term parameters can be a different form, so don't do validation in this case
      #TODO: Change this to do proper validation.
      #if(length(prior_params[[2]][[1]]) != d || length(prior_params[[2]][[2]]) != d){
      if((length(prior_params[[2]][[1]]) != d || length(prior_params[[2]][[2]]) != d) && nm_var != "ind_st_params"){
        stop(paste0("Invalid inverse-gamma parameters for ", nm, " priors (the second element of the ", nm_var, " parameter). This should be a list of length 2, with each entry a numeric of length 1 or ", d, " specifying the shape and scale parameters of the inverse- gamma distributions."))
      }

    }
  return(prior_params)
}

