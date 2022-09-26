#' @rdname PriorConstructorFunctions
#' @export
IndLTPrior <- function(parametrisation_form = "lkj", var_params = list(1, 1), cor_params = 1){
  # validate_parametrisation_form(parametrisation_form)
  # if(is.numeric(cor_params))
  #   cor_params <- list(cor_params)

  ret <- new('IndLTPrior',
             parametrisation_form = parametrisation_form,
             var_params = var_params,
             cor_params = cor_params)
return(ret)
}

setClassUnion(name = "listOrNumeric", members=c("list", "numeric"))

#### Class definition ####
#' A class to hold the priors for the ensemble model.
#'
#' An `IndLTPrior` object encapsulates the prior information for the long-term discrepancies of the individual simulators within the ensemble model.
#' @slot parametrisation_form The parametrisation by which the covariance matrix of the noise of the AR process is decomposed.
#' @slot var_params The parameters characterising the variance of the distribution of the individual long-term discrepancy.
#' @slot cor_params The parameters characterising the correlations of the distribution of the individual long-term discrepancy.
#' @seealso \code{\linkS4class{IndSTPrior}}, \code{\linkS4class{ShaSTPrior}}, \code{\linkS4class{TruthPrior}},
#' @details
#' Individual long-term discrepancies \eqn{\gamma_k} are drawn from a distribtion \deqn{\gamma_k \sim N(0, C_\gamma).} Accepted parametrisation forms for this discrepancy are `lkj`, `beta`, or `inv_wishart`. See details of the `EnsemblePrior()` constructor for more details.
#' @export
setClass(
  "IndLTPrior",
  slots = c(parametrisation_form = "character",
            var_params = "list",
            cor_params = "listOrNumeric")
)



generate_priors_stan_input_ind_lt <- function(d, x){

  if(!is(x, "IndLTPrior")){
    stop("Invalid object for individual long-term priors. This should be a IndLTPrior object.")
  }

  var_priors <- repeat_variance_priors(d, x)

  cor_form   <- correlation_form_prior(x@parametrisation_form)
  cor_params <- correlation_prior(x@cor_params, cor_form, "prior_ind_lt_cor", FALSE)

  ret <- append(list(
      prior_ind_lt_var_a = var_priors[[1]],
      prior_ind_lt_var_b = var_priors[[2]],
      form_prior_ind_lt = cor_form),
    cor_params)

  return(ret)
}

