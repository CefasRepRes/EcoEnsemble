#'@param d A `numeric` giving the number of variables of interest in the ensemble model.
#'@param initial_mean A `numeric` giving the mean of the normal distribution giving the prior on the initial value of the random walk. This is the same value for each variable Default value is `0`.
#'@param initial_var A `numeric` giving the variance of the normal distribution giving the prior on the initial value of the random walk. This is the same value for each variable Default value is `100`.
#'@param rw_covariance A `list` of length `2` containing the inverse-Wishart parameters for the covariance of the random walk of the truth. The default value is `list(2*d, diag(d))`.
#' @examples
#'
#' ##### TruthPrior
#' #Simple default truth prior with 7 variables of interest
#' truth_def <- TruthPrior(7)
#' # A more fine-tuned truth prior for an ensemble with 7 species.
#' truth_cus <- TruthPrior(7, initial_mean = 2, initial_var = 10, rw_covariance = list(10, diag(7)))
#'
#' @rdname PriorConstructorFunctions
#' @export
TruthPrior <- function(d, initial_mean = 0, initial_var = 100, rw_covariance = list(2*d, diag(d))){

  ret <- new('TruthPrior',
             d = d,
             initial_mean  = initial_mean,
             initial_var  = initial_var,
             rw_covariance = rw_covariance)
  return(ret)
}



#### Class definition ####
#' A class to hold the priors for the truth model in the ensemble framework
#'
#' An `TruthPrior` object encapsulates the prior information for the short-term discrepancies of the shared discrepancy of the ensemble model.
#' @slot d A `numeric` giving the number of variables of interest in the ensemble model.
#' @slot initial_mean A `numeric` giving the standard deviation of the normal prior on the initial mean value of the random walk. This is the same standard deviation for each variable of interest.
#' @slot initial_var A `list` of length `2` containing the shape and scale parameters (respectively) for the gamma priors on the variance of the initial value of the truth.
#' @slot rw_covariance A `list` of length `2` containing the inverse-Wishart parameters for the covariance of the random walk of the truth.
#' @details
#' The truth \eqn{\mathbf{y}^{(t)}} is modelled as a random walk such that \deqn{\mathbf{y}^{(t+1)} \sim N(\mathbf{y}^{(t)}, \Lambda_y).} The covariance matrix \eqn{\Lambda_y} is parameterised by an inverse Wishart distribution (contained in the `rw_covariance` slot) and the initial value is modelled as drawn from a normal distribution.
#' @export
setClass(
  "TruthPrior",
  slots = c(d = "numeric",
            initial_mean  = "numeric",
            initial_var  = "numeric",
            rw_covariance = "list")
)



generate_priors_stan_input_truth <- function(d, x){
  if(!is(x, "TruthPrior")){
    stop("Invalid object for priors on the truth. This should be a TruthPrior object.")
  }


  if(length(x@initial_mean) == 1)
    prior_y_init_mean <- rep(x@initial_mean, d)
  if(length(x@initial_var) == 1){
    prior_y_init_var = rep(x@initial_var, d)
  }

  ret <- list(
    prior_y_init_mean = prior_y_init_mean,
    prior_y_init_var = prior_y_init_var,
    prior_sigma_t_inv_wish_nu = x@rw_covariance[[1]],
    prior_sigma_t_inv_wish_sigma = x@rw_covariance[[2]]
  )
  return(ret)
}
