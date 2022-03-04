# These are the currently accepted priors for correlation matrices. If you add
# anything here, you will also need to change the validate_correlation_priors
# function in validation_functions.R
CORRELATIONS_PRIOR_LKJ <- 0
CORRELATIONS_PRIOR_INV_WISHART <- 1
CORRELATIONS_PRIOR_BETA <- 2

correlation_form_prior <- function(prior_choice){
  return(
    switch(EXPR = tolower(prior_choice),
           "lkj" = CORRELATIONS_PRIOR_LKJ,
           "inv_wishart" = CORRELATIONS_PRIOR_INV_WISHART,
           "beta" = CORRELATIONS_PRIOR_BETA)
  )
}

correlation_prior <- function(params, form, type){
  retv <- rep(NA, 5)
  names(retv) <- c(paste0(type, "_lkj"), paste0(type, "_wish_nu"),
                   paste0(type, "_wish_sigma"),paste0(type, "_beta_1"),
                   paste0(type, "_beta_2"))
  ret <- as.list(retv)

  ret[[paste0(type, "_lkj")]] <- numeric(0)
  ret[[paste0(type, "_wish_nu")]] <- numeric(0)
  ret[[paste0(type, "_wish_sigma")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_1")]] <- matrix(0,0,0)
  ret[[paste0(type, "_beta_2")]] <- matrix(0,0,0)


  if (form == CORRELATIONS_PRIOR_LKJ){
    ret[[paste0(type, "_lkj")]] <- array(params, 1)
  }else if(form == CORRELATIONS_PRIOR_INV_WISHART){
    ret[[paste0(type, "_wish_nu")]] <-array(params[[1]], 1)
    ret[[paste0(type, "_wish_sigma")]] <- params[[2]]
  }else if(form == CORRELATIONS_PRIOR_BETA){
    ret[[paste0(type, "_beta_1")]] <- params[[1]]
    ret[[paste0(type, "_beta_2")]] <- params[[2]]
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



# If only one of each inverse-gamma variable is provided, then repeat it for every species considered. We don't do any validation here
repeat_priors <- function(d, prior_params){
    if (is.list(prior_params) && length(prior_params) == 3 && is.list(prior_params[[2]]) && length(prior_params[[2]]) == 2 && length(prior_params[[2]][[1]]) == 1 && length(prior_params[[2]][[2]]) == 1){
      prior_params[[2]][[1]] <- rep(prior_params[[2]][[1]], d)
      prior_params[[2]][[2]] <- rep(prior_params[[2]][[2]], d)
    }
  return(prior_params)
}

