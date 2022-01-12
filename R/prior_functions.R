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
  form_prior_ind_st_cor <- correlation_form_prior(prior_correlations$ind_st_cor_form)
  form_prior_ind_lt_cor <- correlation_form_prior(prior_correlations$ind_lt_cor_form)
  form_prior_sha_st_cor <- correlation_form_prior(prior_correlations$sha_st_cor_form)

  priors_data <- list(
    form_prior_ind_st_cor = form_prior_ind_st_cor,
    form_prior_ind_lt_cor = form_prior_ind_lt_cor,
    form_prior_sha_st_cor = form_prior_sha_st_cor
  )

  priors_data <- append(priors_data, correlation_prior(prior_correlations$ind_st_cor_params, form_prior_ind_st_cor, "prior_ind_st_cor"))
  priors_data <- append(priors_data, correlation_prior(prior_correlations$ind_lt_cor_params, form_prior_ind_lt_cor, "prior_ind_lt_cor"))
  priors_data <- append(priors_data, correlation_prior(prior_correlations$sha_st_cor_params, form_prior_sha_st_cor, "prior_sha_st_cor"))

  return(priors_data)
}

#' Ensemble Discrepancy Priors
#'
#' Encodes prior information on discrepancies for use in the ensemble model.
#' See the supporting information for Spence et. al. (2018) for a worked example of this.
#'
#' @param ind_st_var_params A `numeric` of length 2 giving the shape
#' and scale parameters of the inverse gamma prior on the variance of the individual
#' short-term discrepancies.
#' @param ind_st_cor_form A `character` describing the form of the correlation prior
#' used on the individual short-term discrepancies. Currently supported forms are
#'  `'lkj'`, `'inv_wishart'`, or `'beta'`. See details below.
#' @param ind_st_cor_params A `list` giving the parameters for the parameters for the prior on the correlation matrix of
#' the individual short-term discrepancies. See details below.
#' @param ind_lt_var_params A `vector` of length 2 giving the shape
#' and scale parameters of the inverse gamma prior on the variance of the individual
#' long-term discrepancies.
#' @param ind_lt_cor_form A `character` describing the form of the correlation prior
#' used on the individual long-term discrepancies. Currently supported forms are
#'  `'lkj'`, `'inv_wishart'`, or `'beta'`. See details below.
#' @param ind_lt_cor_params A `list` giving the parameters for the parameters for the prior on the correlation matrix of
#' the individual long-term discrepancies. See details below.
#' @param sha_st_var_exp The exponential parameter for the prior on the variance of
#' the shared short-term discrepancies.
#' @param sha_st_cor_form A `character` describing the form of the correlation prior
#' used on the shared short-term discrepancies. Currently supported forms are
#'  `'lkj'`, `'inv_wishart'`, or `'beta'`. See details below.
#' @param sha_st_cor_params A `list` giving the parameters for the prior on the correlation matrix of
#' the shared short-term discrepancies. See details below.
#' @param sha_lt_sd A `numeric` of the standard deviation of the shared long-term discrepancy. Each
#' element corresponds to the standard deviation of each feature.
#'
#' @details
#' As in Spence et. al. (2018) , the discrepancy covariance matrices (individual and shared short-term discrepancies
#' as well as individual long-term discrepancies) are decomposed
#' into
#' \deqn{\Sigma = \sqrt{\mathrm{diag}(\pi_i)} \Lambda \sqrt{\mathrm{diag}(\pi)},}
#' where \eqn{\pi} is the vector of variances of each species, and \eqn{\Lambda}
#' is the correlation matrix. The variance terms \eqn{\pi}
#' are parameterised by inverse-gamma distributions, and passed through as `ind_st_var_params`,
#' `ind_lt_var_params`, `sha_st_var_params` parameters.
#' @section Correlation matrix priors:
#' There are currently 3 supported prior distributions on correlation matrices which are chosen
#' via the `ind_st_cor_form`, `ind_lt_cor_form`, and `sha_st_cor_form` options.
#' These can take the value `'lkj'`, `'inv_wishart'`, or `'beta'` to
#' refer to the LKJ, inverse Wishart, or Beta distributions respectively. In each case,
#' the associated parameters should be passed through using the relevant `..._cor_params` variable.
#'  * LKJ - The parameter should be a single scalar value \eqn{\eta} giving the LKJ shape parameter as described in the
#'  \href{https://mc-stan.org/docs/2_28/functions-reference/lkj-correlation.html}{Stan manual.} The
#'  density of a correlation matrix is given by \deqn{f(\Sigma | \eta) \alpha det (\Sigma)^{\eta - 1}}
#'  * Inverse Wishart - A `list` containing a scalar value \eqn{\nu} (giving the degrees of
#'  freedom) and a symmetric, positive definite matrix \eqn{\Sigma} (giving the scale
#'  matrix). See the
#'  \href{https://mc-stan.org/docs/2_28/functions-reference/inverse-wishart-distribution.html}{Stan manual}
#'  for more information. The dimensions of \eqn{\Sigma} should be the same as the correlation matrix
#'  it produces (i.e \eqn{N \times N} where \eqn{N} is the number of species)
#'  * Beta - A `list` containing two \eqn{N \times N} matrices (where \eqn{N} is the number of species), giving the prior success parameters \eqn{\alpha}
#'  and prior failure parameters \eqn{\beta} respectively. To ensure positive-definiteness, the
#'  correlations are rescaled from \eqn{[-1,1] \rightarrow [0,1]} via the function
#'  \eqn{\frac{1}{\pi} \tan^{-1} \frac{\rho}{\sqrt{1-\rho^2} + 1/2}}. It is on these
#'  rescaled parameters that the Beta distribution applies.
#'
#' @return A `list` encoding prior information on discrepancies.
#'
#' @references Spence et. al. (2018). A general framework for combining ecosystem models. \emph{Fish and Fisheries}, 19(6):1031-1042.
#' @references Chandler RE. 2013 Exploiting strength, discounting weakness: combining information from multiple climate simulators. \emph{Phil Trans R Soc A} 371: 20120388
#' @examples
#' #Basic usage of function for a model with 4 species.
#' N_species <- 4
#' priors <- define_priors(ind_st_var_params = list(25, 0.25),
#'                         ind_st_cor_form = "lkj", #Using an LKJ distribution for individual short-term discrepancies
#'                         ind_st_cor_params = 30, #The parameter is 30
#'                         ind_lt_var_params = list(rep(25,N_species),rep(0.25,N_species)),
#'                         ind_lt_cor_form = "beta", #For long term discrepancies we use a Beta distribution
#'                         ind_lt_cor_params = list(matrix(40,N_species, N_species), matrix(40, N_species, N_species)),
#'
#'                         sha_st_var_exp = 3,
#'                         sha_st_cor_form = "lkj",
#'                         sha_st_cor_params = 30,
#'                         sha_lt_sd = rep(4,N_species))
#'
#' @export
define_priors <- function(ind_st_var_params, ind_st_cor_form, ind_st_cor_params,
                          ind_lt_var_params, ind_lt_cor_form, ind_lt_cor_params,
                          sha_st_var_exp, sha_st_cor_form, sha_st_cor_params, sha_lt_sd){

  return(
    list(prior_ind_st_var_a = ind_st_var_params[[1]],
       prior_ind_st_var_b = ind_st_var_params[[2]],
       prior_ind_lt_var_a = ind_lt_var_params[[1]],
       prior_ind_lt_var_b = ind_lt_var_params[[2]],
       prior_sha_st_var_exp = sha_st_var_exp,
       prior_sha_lt_sd = sha_lt_sd,
       prior_correlations = list(
         ind_st_cor_form = ind_st_cor_form,
         ind_st_cor_params = ind_st_cor_params,
         ind_lt_cor_form = ind_lt_cor_form,
         ind_lt_cor_params = ind_lt_cor_params,
         sha_st_cor_form = sha_st_cor_form,
         sha_st_cor_params = sha_st_cor_params
         )
       )
    )
  }

