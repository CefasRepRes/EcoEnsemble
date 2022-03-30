#' Constructor for the `EnsemblePrior` class
#'
#' A constructor for the `EnsemblePrior` class. This is used to encode prior information for the ensemble model.
#' @param d A `numeric` specifying the number of species in the ensemble.
#' @param ind_st_params A `list` containing a prior specification for the individual short-term discrepancies \eqn{z_k^{(t)}}. See details
#' @param ind_lt_params A `list` containing a prior specification for the individual long-term discrepancies \eqn{\gamma_k}. See details
#' @param sha_st_params A `list` containing a prior specification for the shared short-term discrepancies \eqn{\eta^{(t)}}. See details
#' @param sha_lt_params A `numeric` containing the standard deviations for the normal prior used on the shared short-term discrepancy \eqn{\mu}. If a single value is supplied,  this is repeated for each species.
#' @param truth_params A `list` containing a prior specification for the processes on the truth \eqn{y^{(t)}}. The default value is `list(10, list(10, 1), list(d, diag(d))` where `d` is the number of species.
#'
#' @details
#' Most discrepancy prior paramaters (`ind_st_params, ind_lt_params, sha_st_params` but not `sha_lt_params`) should be encoded by a `list.` The entries of the list are respectively
#' 1. A `character` specifying how the priors are encoded. Currently supported priors are `'LKJ'`, `'inv_wishart'`, or `'Beta'`.
#' 2. A `list` of length `2` specifying the inverse-gamma parameters for the prior on the discrepancy variances. The first element should be a `numeric` vector containing shape parameters for each species and the second element should be a `numeric` vector containing scale parameters for each species. If only one value is provided in each vector, then this is repeated for all species. See below.
#' 3. A `list` containing the correlation matrix parameters. See below.
#'
#' There are currently three supported prior distributions on covariance matrices. As in Spence et. al. (2018) , the discrepancy covariance matrices (individual and shared short-term discrepancies \eqn{z_k^{(t)},\eta^{(t)}} as well as individual long-term discrepancies \eqn{\gamma_k}, but not the shared long-term discrepancy \eqn{\delta}) are decomposed into a vector of variances and a correlation matrix \deqn{\Lambda = \sqrt{\mathrm{diag}(\pi)}  P \sqrt{\mathrm{diag}(\pi)},} where \eqn{\pi} is the vector of variances for each species, and \eqn{P} is the correlation matrix. The variance terms \eqn{\pi} are parameterised by inverse-gamma distributions, while the correlation matrices can be parameterised in three different ways depending on the choice of parameterisation given in the first element of the `list`. Selecting `'lkj'`, `'inv_wishart'`, or `'beta'` refers to the LKJ, inverse Wishart, or Beta distributions respectively. The associated parameters should be passed through as the third element in the `list`, and should be of the form:
#'  * If `'lkj'` is selected., then  the third element should be a `numeric` \eqn{\eta} giving the LKJ shape parameter, such  that the probability density is given by  (Lewandowski et. al. 2009) \deqn{f(\Sigma | \eta)\propto \mathrm{det} (\Sigma)^{\eta - 1}.}
#'  * If `'inv_wishart'` is selected., then  the third element should be a `list` containing a scalar value \eqn{\nu} (giving the degrees of freedom) and a symmetric, positive definite matrix \eqn{\Sigma} (giving the scale matrix). The dimensions of \eqn{\Sigma} should be the same as the correlation matrix it produces (i.e \eqn{d \times d} where \eqn{d} is the number of species). The density of an inverse Wishart is given by  \deqn{f(W|\eta, S) = \frac{1}{2^{\eta d/2} \Gamma_N \left( \frac{\eta}{2} \right)} |S|^{\eta/2} |W|^{-(\eta + d + 1)/2}  \exp \left(- \frac{1}{2} \mathrm{tr}\left(SW^{-1} \right) \right),} where \eqn{\Gamma_N} is the multivariate gamma function and \eqn{\mathrm{tr \left(X \right)}} is the trace of \eqn{X}.  Note that inverse Wishart distributions act over the space of all covariance matrices. When used for a correlation  matrix, only the subset of valid covariance matrices that are also valid correlation matrices are considered.
#'  * If `'beta'` is selected., then  the third element should be a  `list` containing two \eqn{d \times d} matrices (where \eqn{d} is the number of species), giving the prior success parameters \eqn{\alpha} and prior failure parameters \eqn{\beta} respectively. To ensure positive-definiteness, the correlations are rescaled from \eqn{[-1,1] \rightarrow [0,1]} via the function \eqn{\frac{1}{\pi} \tan^{-1} \frac{\rho}{\sqrt{1-\rho^2}} + \frac{1}{2}}. It is on these rescaled parameters that the Beta distribution applies.
#'
#'
#'  In addition to priors on the discrepancy terms, it is also possible to add prior information on the truth. We require priors on the truth at \eqn{t=0}. By default, a \eqn{N(0, 10)} prior is used on the initial values, and an Inv-Gamma\eqn{(10, 1)} prior is used for the initial variances, however these values can be configured by the `truth_params` argument. The covariance matrix of the random walk of the truth \eqn{\Lambda_y} can also be configured subject to an inverse-Wishart prior. The elements of the `truth_params` list should be
#'  1. A `numeric` giving the standard deviation of the normal prior used for each species of the initial truth values. Default value is 10
#'  2. A `list` of length `2` containing the shape and scale parameters (respectively) for the inverse gamma priors on the initial variance of the truth. The default value is `list(10, 1)`.
#'  3. A `list` of length `2` containing the inverse-Wishart parameters for the random walk of the truth. The default value is `list(d, diag(d))` where `d` is the number of species.
#'
#'
#' @references Spence et. al. (2018). A general framework for combining ecosystem models. \emph{Fish and Fisheries}, 19(6):1031-1042.
#' @examples
#' #Defining priors for a model with 4 species.
#' num_species <- 4
#' priors <- EnsemblePrior(
#'     d = num_species,
#'     ind_st_params = list("lkj",  list(3, 2), 3),
#'     ind_lt_params = list(
#'        "beta",
#'        list(c(10,4,8, 7),c(2,3,1, 4)),
#'        list(matrix(5, num_species, num_species),
#'             matrix(0.5, num_species, num_species))
#'      ),
#'     sha_st_params = list("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'     sha_lt_params = 5,
#'     truth_params = list(10, list(3, 3), list(10, diag(num_species)))
#' )
#' @export
EnsemblePrior <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params = list(10, list(10, 1), list(d, diag(d)))){
  validate_input(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params)
  priors_stan_input <- create_prior_stan_input(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params)
  ret <- new('EnsemblePrior',
             d = d,
             ind_st_params = ind_st_params,
             ind_lt_params = ind_lt_params,
             sha_st_params = sha_st_params,
             sha_lt_params = sha_lt_params,
             truth_params = truth_params,
             priors_stan_input = priors_stan_input)
  return(ret)
}



#### Class definition ####
#' A class to hold the priors for the ensemble model.
#'
#' An `EnsemblePrior` object encapsulates the prior information for the ensemble model.
#'
#' @slot d A `numeric` specifying the number of species in the ensemble.
#' @slot ind_st_params A `list` containing a prior specification for the individual short-term discrepancies \eqn{z_k^{(t)}}. See details
#' @slot ind_lt_params A `list` containing a prior specification for the individual long-term discrepancies \eqn{\gamma_k}. See details
#' @slot sha_st_params A `list` containing a prior specification for the shared short-term discrepancies \eqn{\eta^{(t)}}. See details
#' @slot sha_lt_params A `numeric` containing the standard deviations for the normal prior used on the shared short-term discrepancy \eqn{\mu}. If a single value is supplied,  this is repeated for each species.
#' @slot truth_params A `list` containing a prior specification for the processes on the truth \eqn{y^{(t)}}. See details. The default value is `list(10, list(10, 1), list(10, d))` where `d` is
#' @slot priors_stan_input A `list` containing the prior data in the correct form to fit the model in Stan. This information is automatically generated by the constructor.
#'
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#' @references Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe. 2009. “Generating Random Correlation Matrices Based on Vines and Extended Onion Method.” Journal of Multivariate Analysis 100: 1989–2001.
#' @seealso \code{\link{EnsemblePrior}}, \code{\link{EnsembleData}}, \code{\link{fit_ensemble_model}},
#' @export
setClass(
  "EnsemblePrior",
  slots = c(d = "numeric",
            ind_st_params = "list",
            ind_lt_params = "list",
            sha_st_params = "list",
            sha_lt_params = "numeric",
            truth_params = "list",
            priors_stan_input = "list"
  )
)





create_prior_stan_input <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params){

  validate_input(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params)

  # If only one item is specified for discrepancy priors, repeat it to be the same value for each species.
  ind_st_params <- repeat_priors(d, ind_st_params, "individual short-term", "ind_st_params")
  ind_lt_params <- repeat_priors(d, ind_lt_params, "individual long-term", "ind_lt_params")
  sha_st_params <- repeat_priors(d, sha_st_params, "shared short-term", "sha_st_params")
  if(is.numeric(sha_lt_params) && length(sha_lt_params) == 1)
    sha_lt_params <- rep(sha_lt_params, d)

  #Do the same for truth parameters
  if(is.list(truth_params) && length(truth_params) == 3){
    y_init_mean <- truth_params[[1]]
    if(length(y_init_mean) == 1)
      truth_params[[1]] <- rep(y_init_mean, d)

    y_init_inv_gamma <- truth_params[[2]]
    if(is.list(y_init_inv_gamma) && length(y_init_inv_gamma) == 2 && length(y_init_inv_gamma[[1]]) == 1 && length(y_init_inv_gamma[[2]]) == 1){
      truth_params[[2]][[1]] <- rep(y_init_inv_gamma[[1]], d)
      truth_params[[2]][[2]] <- rep(y_init_inv_gamma[[2]], d)
    }

  }

  priors <- list(
    prior_ind_st_var_a = ind_st_params[[2]][[1]],
    prior_ind_st_var_b = ind_st_params[[2]][[2]],
    prior_ind_lt_var_a = ind_lt_params[[2]][[1]],
    prior_ind_lt_var_b = ind_lt_params[[2]][[2]],
    prior_sha_st_var_a = sha_st_params[[2]][[1]],
    prior_sha_st_var_b = sha_st_params[[2]][[2]],
    prior_sha_lt_sd = sha_lt_params,
    prior_y_init_mean_sd = truth_params[[1]],
    prior_y_init_var_a = truth_params[[2]][[1]],
    prior_y_init_var_b = truth_params[[2]][[2]],
    prior_sigma_t_inv_wish_nu = truth_params[[3]][[1]],
    prior_sigma_t_inv_wish_sigma = truth_params[[3]][[2]]
  )
  priors_correlations <- list(
    ind_st_cor_form = ind_st_params[[1]],
    ind_st_cor_params = ind_st_params[[3]],
    ind_lt_cor_form = ind_lt_params[[1]],
    ind_lt_cor_params = ind_lt_params[[3]],
    sha_st_cor_form = sha_st_params[[1]],
    sha_st_cor_params = sha_st_params[[3]]
  )
  priors <- append(priors, generate_correlation_priors_stan_data(priors_correlations))
  return(priors)
}



validate_input <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params){
  if(!is.numeric(d) || length(d) != 1)
    stop(paste0("The number of species (d) should be an integer value, current value: ", d))

  validate_priors_ar_process(ind_st_params, "individual short-term", d)
  validate_priors_ar_process(ind_lt_params, "individual long-term", d)
  validate_priors_ar_process(sha_st_params, "shared short-term", d)

}
