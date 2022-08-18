validate_prior_compatibility <- function(type_str, params, d) {

  var_params <- params@var_params
  cor_params <- params@cor_params
  form <- tolower(params@parametrisation_form)

  msg <- ""
  if(form == "lkj" &
     length(cor_params) != 1){
    stop("Invalid lkj parameter for ", type_str, " correlation matrix prior. This should be a list with exactly one real number.")
  }else if (form == "inv_wishart"){
    if(length(cor_params) != 2 || !is.numeric(cor_params[[1]]) || !is.matrix(cor_params[[2]]) || !is.square.matrix(cor_params[[2]]) || cor_params[[1]] <= d-1){
      stop("Invalid inverse Wishart parameters for ", type_str, " correlation matrix prior. This should be a numeric giving the degrees of freedom and a d x d scale matrix, with degrees of freedom > d-1.")
    }
  }else if (form == "beta" &&
            (length(cor_params) !=2 ||
             !is.matrix(cor_params[[1]]) ||
             !is.matrix(cor_params[[2]]) ||
             !is.square.matrix(cor_params[[1]]) ||
             !is.square.matrix(cor_params[[2]]) ||
             !is.symmetric.matrix(cor_params[[1]]) ||
             !is.symmetric.matrix(cor_params[[2]]) ||
             !(dim(cor_params[[1]])[1] == d) ||
             !(dim(cor_params[[2]])[1] == d))){
    stop("Invalid beta parameters for ", type_str, " correlation matrix priors. These should be square, symmetric matrices with the same dimension as the number of model outputs.")
  }else if (form == "hierarchical"& length(cor_params) !=4){
    stop("Invalid parameters for ", type_str, " correlation matrix priors. This should be a list of length 4.")
  }
}



validate_parametrisation_form <- function(form, valid_forms = c(CORRELATIONS_PRIOR_LKJ,CORRELATIONS_PRIOR_INV_WISHART,CORRELATIONS_PRIOR_BETA)){
  if(!isTRUE(correlation_form_prior(form) %in% valid_forms)){
    stop("Invalid parametrisation choice for priors. Prior parametrisation forms should be one of 'lkj', 'inv_wishart', 'beta', or 'hierarchical'. The 'hierarchical' option is only available for individual short-term discrepancies. Prior choice: ", form)
  }

}


validate_prior_AR_params <- function(AR_params){
  if(!is.numeric(AR_params) || !(length(AR_params) == 2)){
    stop("Invalid autoregressive parameters specified. This should be a numeric of length 2 giving the shape parameters of the Beta distribution used as a prior for the AR parameters.")
  }
}
