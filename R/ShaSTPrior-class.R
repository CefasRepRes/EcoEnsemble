#' @rdname PriorConstructorFunctions
#' @export
ShaSTPrior <- function(parametrisation_form = "lkj", var_params = list(1, 1), cor_params = 1, AR_params = c(1,1)){

  validate_parametrisation_form(parametrisation_form)
  validate_prior_AR_params(AR_params)

  # if(is.numeric(cor_params))
  #   cor_params <- list(cor_params)


  ret <- new('ShaSTPrior',
             AR_param = AR_params,
             parametrisation_form = parametrisation_form,
             var_params = var_params,
             cor_params = cor_params)
  return(ret)
}



#### Class definition ####
#' A class to hold the priors for the ensemble model.
#'
#' An `ShaSTPrior` object encapsulates the prior information for the short-term discrepancies of the shared discrepancy of the ensemble model.
#' @slot AR_param The parameters giving the Beta parameters for the main parameter of the AR(1) process.
#' @slot parametrisation_form The parametrisation by which the covariance matrix of the noise of the AR process is decomposed.
#' @slot var_params The parameters characterising the variance of the AR process on the shared short-term discrepancy.
#' @slot cor_params The parameters characterising the correlations of the AR process on the shared short-term discrepancy.
#' @details
#' Shared short-term discrepancies \eqn{\mathbf{\eta}^{(t)}} are modelled as an AR(1) process so that \deqn{\mathbf{\eta}^{(t+1)} \sim N(R_\eta \mathbf{\eta}, \Lambda_\eta).} Accepted parametrisation forms for this discrepancy are `lkj`, `beta`, or `inv_wishart`. See details of the `EnsemblePrior()` constructor for more details.
#' @export
setClass(
  "ShaSTPrior",
  slots = c(AR_param = "numeric",
            parametrisation_form = "character",
            var_params = "list",
            cor_params = "listOrNumeric")
)



generate_priors_stan_input_sha_st <- function(d, x){
  if(!is(x, "ShaSTPrior")){
    stop("Invalid object for shared short-term priors. This should be a ShaSTPrior object.")
  }

  var_priors <- repeat_variance_priors(d, x)

  cor_form   <- correlation_form_prior(x@parametrisation_form)
  cor_params <- correlation_prior(x@cor_params, cor_form, "prior_sha_st_cor", FALSE)

  ret <- append(list(
    prior_sha_st_var_a = var_priors[[1]],
    prior_sha_st_var_b = var_priors[[2]],
    prior_sha_st_ar_alpha = x@AR_param[1],
    prior_sha_st_ar_beta = x@AR_param[2],
    form_prior_sha_st = cor_form),
    cor_params)

  return(ret)

}
