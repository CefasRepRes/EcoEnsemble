#' @param parametrisation_form The parametrisation by which the covariance matrix of the noise of the AR process (in the case of `IndSTPrior` and `ShaSTPrior` objects or the covariance of the distribution of long-term discrepancies for `IndLTPrior` objects) is decomposed. The default is `hierarchical` for `IndSTPrior` objects, and `lkj` otherwise. See details.
#' @param var_params The parameters characterising the variance of the AR process (in the case of `IndSTPrior` and `ShaSTPrior` objects or the variance of the distribution of long-term discrepancies for `IndLTPrior` objects) on the discrepancy. The default value is `list(-3, 1, 8, 4)` for `IndSTPrior` objects, `list(1, 1)` for `IndLTPrior` objects, and `list(1, 10)` for `ShaSTPrior` objects. See details.
#' @param cor_params The parameters characterising the correlations of the AR process (or the distribution of long-term discrepancies) on the short-term discrepancies. The default value in this case is to use `list(0.1,0.1,0.1,0.1)` for `IndSTPrior` objects, and `1` for `IndLTPrior` and `ShaSTPrior` objects. See details.
#' @param AR_params The parameters giving the beta parameters for the prior distribution on the autoregressive parameter of the AR(1) process. The default is `c(2,2)`. See details.
#' @examples
#' ##### Different forms of the individual long term discrepancy priors
#' #LKJ(10) priors on correlation matrices and gamma(5, 3) priors on the variances
#' ist_lkj <- IndSTPrior("lkj", list(5, 3), 10)#
#'
#' #Same as above but with an additional beta(2, 4) prior on
#' #the autoregressive parameter of the AR process.
#' ist_lkj <- IndSTPrior("lkj", list(5, 3), 10, AR_params = c(2, 4))
#'
#' #Same as above but with different variance priors for 5 different variables of interest.
#' #This encodes that there is a gamma(1, 1) prior on the variance of the first variable,
#' #a gamma(23, 1) on the second variable etc...
#' ist_lkj <- IndSTPrior("lkj", list(c(1,23,24,6,87), c(1,1,1,1,5)), 10, AR_params = c(2, 4))
#'
#' #Hierarchical priors with gamma(1,2) and gamma(10, 1) on the variance hyperparameters and
#' #gamma(3,4), gamma(5,6) on the correlation hyperparameters
#' ist_hie <- IndSTPrior("hierarchical", list(1,2,10,1), list(3,4,5,6))
#'
#' #Hierarchical priors with gamma(1,2) and gamma(10, 1) on the variance hyperparameters and
#' #the beta conjugate prior with parameters (p = 0.75, q = 0.75, k = 0.2) on the correlation hyperparameters
#' ist_hie_beta_conj <- IndSTPrior("hierarchical_beta_conjugate", list(1,2,10,1), list(0.75, 0.75, 0.2))
#'
#' #Inverse Wishart correlation priors. Gamma(2, 1/3) priors are on the variances and
#' #inv-Wishart(5, diag(5)) on the correlation matrices.
#' ist_inW <- IndSTPrior("inv_wishart", list(2, 1/3),list(5, diag(5)))
#' @rdname PriorConstructorFunctions
#' @export

IndSTPrior <- function(parametrisation_form = "hierarchical", var_params= list(-3, 1, 8, 4), cor_params = list(0.1, 0.1, 0.1, 0.1), AR_params = c(2, 2)){

  validate_parametrisation_form(parametrisation_form, valid_forms = c(CORRELATIONS_PRIOR_LKJ,
                                                                      CORRELATIONS_PRIOR_INV_WISHART,
                                                                      CORRELATIONS_PRIOR_BETA,
                                                                      CORRELATIONS_PRIOR_HIERARCHICAL,
                                                                      CORRELATIONS_PRIOR_BETA_CONJUGATE))
  validate_prior_AR_params(AR_params)

  # if(is.numeric(cor_params))
  #   cor_params <- list(cor_params)

  ret <- new('IndSTPrior',
             AR_param = AR_params,
             parametrisation_form = parametrisation_form,
             var_params = var_params,
             cor_params = cor_params)
  return(ret)
}



#### Class definition ####
#' A class to hold the priors for the ensemble model.
#'
#' An `IndSTPrior` object encapsulates the prior information for the short-term discrepancies of the individual simulators within the ensemble model.
#' @slot AR_param The parameters giving the beta parameters for the autoregressive parameter of the AR(1) process.
#' @slot parametrisation_form The parametrisation by which the covariance matrix of the noise of the AR process is decomposed.
#' @slot var_params The parameters characterising the variance of the AR process on the individual short-term discrepancy.
#' @slot cor_params The parameters characterising the correlations of the AR process on the individual short-term discrepancy. .
#' @seealso \code{\linkS4class{IndLTPrior}}, \code{\linkS4class{ShaSTPrior}}, \code{\linkS4class{TruthPrior}},
#' @details
#' Individual short-term discrepancies \eqn{z_k^{(t)}} are modelled as an AR(1) process so that \deqn{z_k^{(t+1)} \sim N(R_k z_k^{(t)}, \Lambda_k).} Accepted parametrisation forms for this discrepancy are `lkj`, `beta`, `inv_wishart`, `hierarchical` or `hierarchical_beta_conjugate`. See details of the `EnsemblePrior()` constructor for more details.
#' @export
setClass(
  "IndSTPrior",
  slots = c(AR_param = "numeric",
            parametrisation_form = "character",
            var_params = "listOrNumeric",
            cor_params = "listOrNumeric")
)




generate_priors_stan_input_ind_st <- function(d, x){

  if(!is(x, "IndSTPrior")){
    stop("Invalid object for individual short-term priors. This should be a IndSTPrior object.")
  }

  cor_form   <- correlation_form_prior(x@parametrisation_form)
  cor_params <- correlation_prior(x@cor_params, cor_form, "prior_ind_st_cor", TRUE)

  var_priors <- list(
    prior_ind_st_var_a = numeric(0),
    prior_ind_st_var_b = numeric(0),
    prior_ind_st_var_hierarchical_hyperparams = numeric(0)
  )
  if(cor_form == CORRELATIONS_PRIOR_HIERARCHICAL || cor_form == CORRELATIONS_PRIOR_BETA_CONJUGATE){
    var_priors$prior_ind_st_var_hierarchical_hyperparams <- unlist(x@var_params)
  }else{
    repeated_var_params <- repeat_variance_priors(d, x)
    var_priors$prior_ind_st_var_a <- repeated_var_params[[1]]
    var_priors$prior_ind_st_var_b <- repeated_var_params[[2]]
  }





  ret <- append(var_priors, cor_params)
  ret <- append(ret,
                list(form_prior_ind_st = cor_form,
                     prior_ind_st_ar_alpha = x@AR_param[1],
                     prior_ind_st_ar_beta = x@AR_param[2]))
  return(ret)

}
