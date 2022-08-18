#'@param d A `numeric` giving the number of variables of interest in the ensemble model.
#'@param initial_mean_sd A `numeric` giving the standard deviation of the normal prior on the initial mean value of the random walk. This is the same standard deviation for each variable Default value is `10`.
#'@param initial_vars A `list` of length `2` containing the shape and scale parameters (respectively) for the inverse gamma priors on the variance of the initial value of the truth. The default value is `list(10, 1)`.
#'@param rw_covariance A `list` of length `2` containing the inverse-Wishart parameters for the covariance of the random walk of the truth. The default value is `list(d, diag(d))`.
#' @examples
#'
#' ##### TruthPrior
#' #Simple default truth prior with 7 variables of interest
#' truth_def <- TruthPrior(7)
#' # A more fine-tuned truth prior for an ensemble with 7 species.
#' truth_cus <- TruthPrior(7, initial_mean_sd = 2, initial_vars = list(1, 1), rw_covariance = list(10, diag(7)))
#'
#' @rdname PriorConstructorFunctions
#' @export
TruthPrior <- function(d, initial_mean_sd = 10, initial_vars = list(10, 1), rw_covariance = list(d, diag(d))){

  ret <- new('TruthPrior',
             d = d,
             initial_mean_sd  = initial_mean_sd,
             initial_vars  = initial_vars,
             rw_covariance = rw_covariance)
  return(ret)
}



#### Class definition ####
#' A class to hold the priors for the truth model in the ensemble framework
#'
#' An `TruthPrior` object encapsulates the prior information for the short-term discrepancies of the shared discrepancy of the ensemble model.
#' @slot d A `numeric` giving the number of variables of interest in the ensemble model.
#' @slot initial_mean_sd A `numeric` giving the standard deviation of the normal prior on the initial mean value of the random walk. This is the same standard deviation for each variable of interest.
#' @slot initial_vars A `list` of length `2` containing the shape and scale parameters (respectively) for the inverse gamma priors on the variance of the initial value of the truth.
#' @slot rw_covariance A `list` of length `2` containing the inverse-Wishart parameters for the covariance of the random walk of the truth.
#' @details
#' The truth \eqn{\mathbf{y}^{(t)}} is modelled as a random walk such that \deqn{\mathbf{y}^{(t+1)} \sim N(\mathbf{y}^{(t)}, \Lambda_y).} The covariance matrix \eqn{\Lambda_y} is parameterised by an inverse Wishart distribution (contained in the `rw_covariance` slot) and the initial value and covariance are modelled as drawn from a normal distribution (for the value) and an inverse gamma distribution (for the variance).
#' @export
setClass(
  "TruthPrior",
  slots = c(d = "numeric",
            initial_mean_sd  = "numeric",
            initial_vars  = "list",
            rw_covariance = "list")
)



generate_priors_stan_input_truth <- function(d, x){
  if(class(x) != "TruthPrior"){
    stop("Invalid object for priors on the truth. This should be a TruthPrior object.")
  }


  if(length(x@initial_mean_sd) == 1)
    prior_y_init_mean_sd <- rep(x@initial_mean_sd, d)
  if(length(x@initial_vars) == 2 && length(x@initial_vars[[1]]) == 1 && length(x@initial_vars[[2]]) == 1){
    prior_y_init_var_a = rep(x@initial_vars[[1]], d)
    prior_y_init_var_b = rep(x@initial_vars[[2]], d)
  }

  ret <- list(
    prior_y_init_mean_sd = prior_y_init_mean_sd,
    prior_y_init_var_a = prior_y_init_var_a ,
    prior_y_init_var_b = prior_y_init_var_b,
    prior_sigma_t_inv_wish_nu = x@rw_covariance[[1]],
    prior_sigma_t_inv_wish_sigma = x@rw_covariance[[2]]
  )
  return(ret)
}
