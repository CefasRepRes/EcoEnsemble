#' Constructor for the `EnsemblePrior` class
#'
#' Constructors for the `EnsemblePrior` class and related classes. These functions are used to encode prior information for the ensemble model. The `IndSTPrior`, `IndLTPrior`, `ShaSTPrior`, and `TruthPrior` constructors encapsulate prior information about the discrepancies.
#' @param d A `numeric` specifying the number of variables of interest in the ensemble.
#' @param ind_st_params An `IndSTPrior` object specifying priors for the individual short-term discrepancies \eqn{z_k^{(t)}}.
#' @param ind_lt_params An `IndLTPrior` object specifying priors for the individual long-term discrepancies \eqn{\gamma_k}.
#' @param sha_st_params A `ShaSTPrior` object specifying priors for the shared short-term discrepancies \eqn{\eta^{(t)}}.
#' @param sha_lt_params A `numeric` of length `d` or `1` containing the standard deviations for the normal prior used on the shared short-term discrepancy \eqn{\mu}. If a single value is supplied,  this is repeated for each variable of interest.
#' @param truth_params A `TruthPrior` object specifying priors for the processes on the truth \eqn{y^{(t)}}. The default value is `TruthPrior(d)`.
#' @details
#'`IndSTPrior` and `ShaSTPrior` discrepancy prior parameter objects contain 4 slots corresponding to:
#' 1. `parametrisation_form` - A `character` specifying how the priors are parametrised. Currently supported priors are `'LKJ'`, `'inv_wishart'`, `'beta'`, or `'hierarchical'` (`'hierarchical'` is only supported for `IndSTPrior` objects).
#' 2. `var_params` - The prior parameters for the discrepancy variances, either a `list` of length `2` or a `numeric` of length `4`. See below.
#' 3. `cor_params` - The correlation matrix parameters, either a `list` of length `2` or a `numeric` of length `4`. See below.
#' 4. `AR_params` - Parameters for the autoregressive parameter as a `numeric` of length `2`.
#'
#' `IndLTPrior` discrepancy prior parameter objects contain the slots `parametrisation_form`, `var_params`, and `cor_params`.
#'
#' There are currently four supported prior distributions on covariance matrices. As in Spence et. al. (2018), the individual and shared short-term discrepancy covariances, \eqn{\Lambda_k} and \eqn{\Lambda_\eta}, as well as the individual long-term discrepancy covariance, \eqn{\Lambda_\gamma},  are decomposed into a vector of variances and a correlation matrix \deqn{\Lambda = \sqrt{\mathrm{diag}(\pi)}  P \sqrt{\mathrm{diag}(\pi)},} where \eqn{\pi} is the vector of variances for each variable of interest (VoI), and \eqn{P} is the correlation matrix.
#'
#'   Selecting `'lkj'`, `'inv_wishart'`, `'beta'`, or `'hierarchical'` refers to setting LKJ, inverse Wishart, beta, or hierarchical prior distributions on the covariance matrix respectively. The variance parameters should be passed through as the `var_params` slot of the object and the correlation parameters should be passed through as the `cor_params`. For `'lkj'`, `'inv_wishart'`, and `'beta'` selections, variances are parameterised by inverse-gamma distributions, so the `var_params` slot should be a `list` of length two, where each element gives the shape and scale parameters for each VoI (either as a single value which is the same for each VoI or a `numeric` with the same length as the number of VoI). The correlations should be in the following form:
#'  * If `'lkj'` is selected, then `cor_params` should be a `numeric` \eqn{\eta} giving the LKJ shape parameter, such  that the probability density is given by  (Lewandowski et. al. 2009) \deqn{f(\Sigma | \eta)\propto \mathrm{det} (\Sigma)^{\eta - 1}.} Variances are parameterised by inverse-gamma distributions.
#'  * If `'inv_wishart'` is selected, then  `cor_params` should be a `list` containing a scalar value \eqn{\nu} (giving the degrees of freedom) and a symmetric, positive definite matrix \eqn{\Sigma} (giving the scale matrix). The dimensions of \eqn{\Sigma} should be the same as the correlation matrix it produces (i.e \eqn{d \times d} where \eqn{d} is the number of VoI). The density of an inverse Wishart is given by  \deqn{f(W|\eta, S) = \frac{1}{2^{\eta d/2} \Gamma_N \left( \frac{\eta}{2} \right)} |S|^{\eta/2} |W|^{-(\eta + d + 1)/2}  \exp \left(- \frac{1}{2} \mathrm{tr}\left(SW^{-1} \right) \right),} where \eqn{\Gamma_N} is the multivariate gamma function and \eqn{\mathrm{tr \left(X \right)}} is the trace of \eqn{X}.  Note that inverse Wishart distributions act over the space of all covariance matrices. When used for a correlation  matrix, only the subset of valid covariance matrices that are also valid correlation matrices are considered. Variances are parameterised by inverse-gamma distributions.
#'  * If `'beta'` is selected, then  `cor_params` should be a  `list` containing two symmetric `d`\eqn{\times}`d` matrices \eqn{A} and \eqn{B} giving the prior success parameters and prior failure parameters respectively. The correlation between the `i`th and `j`th VoI is \eqn{\rho_{i, j}} with \deqn{\frac{1}{\pi} \tan^{-1} \frac{\rho_{i, j}}{\sqrt{1-\rho_{i, j}^2}} + \frac{1}{2} \sim \mathrm{beta}(A_{i, j}, B_{i, j}).} Variances are parameterised by inverse-gamma distributions.
#'  * If `'hierarchical'` is selected, then variances are parameterised by hierarchical gamma distributions:
#'  \deqn{\pi_{k, i} \sim \mathrm{gamma}(a_{k, i}, b_{k, i})} with priors
#'  \deqn{a_{k, i} \sim \mathrm{gamma}(\alpha_\pi, \beta_\pi),}
#'  \deqn{b_{k, i} \sim \mathrm{gamma}(\gamma_\pi, \delta_\pi).}
#'  The `var_params` slot should then be a `numeric` of length 4, giving the \eqn{\alpha_\pi, \beta_\pi, \gamma_\pi, \delta_\pi} hyperparameters respectively. Correlations (\eqn{\rho_{k, i, j}} where \eqn{\rho_{k, i, j}} is the correlation between VoI \eqn{i} and \eqn{j} for the \eqn{k}th simulator) are parameterised by hierarchical beta distributions.
#'  \deqn{\frac{\rho_{k, i, j} + 1}{2} \sim \mathrm{beta}(c_{k, i, j}, d_{k, i, j})} with priors
#'  \deqn{c_{k, i, j} \sim \mathrm{gamma}(\alpha_\rho, \beta_\rho),}
#'  \deqn{d_{k, i, j} \sim \mathrm{gamma}(\gamma_\rho, \delta_\rho).}
#'  The `cor_params` slot should be a `numeric` of length 4 giving the \eqn{\alpha_\rho, \beta_\rho, \gamma_\rho, \delta_\rho} hyperparameters. respectively. NOTE: This options is only supported for the individual short-term discrepancy terms.
#'
#'Priors may also be specified for the autoregressive parameters for discrepancies modelled using autoregressive processes (i.e. for `IndSTPrior` and `ShaSTPrior` objects). These are parametrised via beta distributions such that the autoregressive parameter \eqn{R \in (-1,1)} satisfies \deqn{\frac{R+1}{2} \sim \mathrm{Beta}(\alpha, \beta)}.
#'
#'  In addition to priors on the discrepancy terms, it is also possible to add prior information on the truth. We require priors on the truth at \eqn{t=0}. By default, a \eqn{N(0, 10)} prior is used on the initial values, and an Inv-Gamma\eqn{(10, 1)} prior is used for the initial variances, however these values can be configured by the `truth_params` argument. The covariance matrix of the random walk of the truth \eqn{\Lambda_y} can also be configured subject to an inverse-Wishart prior. The `truth_params` argument should be a `TruthPrior` object.
#' @references Spence et. al. (2018). A general framework for combining ecosystem models. \emph{Fish and Fisheries}, 19(6):1031-1042.
#'
#' @examples
#' #Defining priors for a model with 4 species.
#' num_species <- 4
#' priors <- EnsemblePrior(
#'   d = num_species,
#'   ind_st_params = IndSTPrior("lkj",  list(3, 2), 3, AR_params = c(1,1)),
#'   ind_lt_params = IndLTPrior(
#'     "beta",
#'     list(c(10,4,8, 7),c(2,3,1, 4)),
#'     list(matrix(5, num_species, num_species),
#'          matrix(0.5, num_species, num_species))
#'   ),
#'   sha_st_params = ShaSTPrior("inv_wishart",list(2, 1/3),list(5, diag(num_species))),
#'   sha_lt_params = 5,
#'   truth_params = TruthPrior(num_species, 10, list(3, 3), list(10, diag(num_species)))
#' )
#'
#' @return `EnsemblePrior` returns an object of class `EnsemblePrior`.
#' `IndSTPrior` returns an object of class `IndSTPrior`.
#' `IndLTPrior` returns an object of class `IndLTPrior`.
#' `ShaSTPrior` returns an object of class `ShaSTPrior`.
#' `TruthPrior` returns an object of class `TruthPrior`.
#' @rdname PriorConstructorFunctions
#' @export
EnsemblePrior <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params = TruthPrior(d)){

  validate_prior(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params)

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
#' @slot d A `numeric` specifying the number of variables of interest in the ensemble.
#' @slot ind_st_params A `list` containing a prior specification for the individual short-term discrepancies \eqn{z_k^{(t)}}. See details of the `EnsemblePrior()` constructor.
#' @slot ind_lt_params A `list` containing a prior specification for the individual long-term discrepancies \eqn{\gamma_k}. See details of the `EnsemblePrior()` constructor.
#' @slot sha_st_params A `list` containing a prior specification for the shared short-term discrepancies \eqn{\eta^{(t)}}. See details of the `EnsemblePrior()` constructor.
#' @slot sha_lt_params A `numeric` containing the standard deviations for the normal prior used on the shared short-term discrepancy \eqn{\mu}. If a single value is supplied,  this is repeated for each variable
#' @slot truth_params A `list` containing a prior specification for the processes on the truth \eqn{y^{(t)}}. See details of the `EnsemblePrior()` constructor. The default value is `TruthPrior(d)`.
#' @slot priors_stan_input A `list` containing the prior data in the correct form to fit the model in Stan. This information is automatically generated by the constructor.
#'
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#' @references Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe. 2009. “Generating Random Correlation Matrices Based on Vines and Extended Onion Method.” Journal of Multivariate Analysis 100: 1989–2001.
#' @export
setClass(
  "EnsemblePrior",
  slots = c(d = "numeric",
            ind_st_params = "IndSTPrior",
            ind_lt_params = "IndLTPrior",
            sha_st_params = "ShaSTPrior",
            sha_lt_params = "numeric",
            truth_params = "TruthPrior",
            priors_stan_input = "list"
  )
)



create_prior_stan_input <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params){

  if(is.numeric(sha_lt_params) && length(sha_lt_params) == 1)
    sha_lt_params <- rep(sha_lt_params, d)

  dat <- generate_priors_stan_input_ind_st(d, ind_st_params) %>%
    append(generate_priors_stan_input_ind_lt(d, ind_lt_params)) %>%
    append(generate_priors_stan_input_sha_st(d, sha_st_params)) %>%
    append(generate_priors_stan_input_truth(d, truth_params)) %>%
    append(list(prior_sha_lt_sd = sha_lt_params))
  return(dat)
}



validate_prior <- function(d, ind_st_params, ind_lt_params, sha_st_params, sha_lt_params, truth_params){
  if(!is.numeric(d) || length(d) != 1)
    stop(paste0("The number of variables of interest (d) should be an integer value, current value: ", d))


  if(!is.numeric(sha_lt_params) || !(length(sha_lt_params) %in% c(1, d)))
    stop("The shared long-term discrepancy parameters should be a numeric of length 1 or d.")


  validate_prior_compatibility("individual short-term",ind_st_params,  d)
  validate_prior_compatibility("individual long-term", ind_lt_params, d)
  validate_prior_compatibility("shared short-term", sha_st_params, d)

}
