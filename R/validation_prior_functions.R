validate_priors_ar_process <- function(ar_priors, nm, d){
  #TODO: Add unit tests for these validation steps.
  if(!is.list(ar_priors) || length(ar_priors) != 3 || !is.character(ar_priors[[1]])){
    msg <- paste0("Invalid prior specification for the ", nm, " discrepancies. This should be a list of length 3 with entries respectively specifying (1)The choice of parametrisation (2)The inverse-gamma parameters for the variance terms and (3) The correlation matrix parametrers. ")
    stop(msg)
  }

  validate_correlation_priors(nm, ar_priors[[1]], ar_priors[[3]], d)
}


validate_correlation_priors <- function(type_str, form, params, d) {
  form <- tolower(form)
  if (!(form %in% c("lkj", "inv_wishart", "beta"))){
    msg <- paste0("Invalid prior choice on ", type_str," correlation matrix. Correlation priors should be one of 'lkj', 'inv_wishart', or 'beta'. Prior choice: ", form)
    stop(msg)
  }

  msg <- ""
  if(form == "lkj" &
     length(params) != 1){
    msg <- paste0("Invalid lkj parameter for ", type_str, " correlation matrix prior. This should be a list with exactly one real number.")

  }else if (form == "inv_wishart"){
    if(length(params) != 2 || !is.numeric(params[[1]]) || !is.matrix(params[[2]]) || !is.square.matrix(params[[2]]) || params[[1]] <= d-1){
      msg <- paste0("Invalid inverse Wishart parameters for ", type_str, " correlation matrix prior. This should be a numeric giving the degrees of freedom and a d x d scale matrix, with degrees of freedom > d-1.")
    }
  }else if (form == "beta" &
            (length(params) !=2 ||
             !is.matrix(params[[1]]) ||
             !is.matrix(params[[2]]) ||
             !is.square.matrix(params[[1]]) ||
             !is.square.matrix(params[[2]]) ||
             !is.symmetric.matrix(params[[1]]) || #TODO:What are the exact conditions we should check here? Do we need symmetry?
             !is.symmetric.matrix(params[[2]]) ||
             !(dim(params[[1]])[1] == d) ||
             !(dim(params[[2]])[1] == d))){
    msg <- paste0("Invalid beta parameters for ", type_str, " correlation matrix priors. These should be square, symmetric matrices with the same dimension as the number of model outputs.")
  }
  if(msg != ""){
    stop(msg)
  }

}
