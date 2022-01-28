#' Constructor for the `EnsembleData` class
#'
#' A constructor for the `EnsembleData` class. This is used to convert
#' input data into the required form to fit the ensemble model using Stan.
#'@param observations A `list` of observations and associated covariances.
#'See details.
#'@param simulators A `list` of simulator output and covariance pairs. See details.
#'@param priors A `list` of data specifying the prior distributions for the ensemble.
#'See the `define_priors` function for the required form.
#'@details Observation and covariance pairs should be passed through as a `list`,
#'with observations / model outputs as the first element, and the covariance matrix
#'as the second. The observations / model outputs should be a `data.frame` (or a `matrix`) with each
#'column giving a different species and each row a different time. Rows should be named with
#'the times and columns should be named the appropriate variables. It is fine to have missing
#'years or species, however for any model output there should be at least one observation in the observation
#'data frame. Missing species / years are encoded by simply not including the corresponding row / column in the
#'data frame / matrix.
#'@return An object of type \linkS4class{EnsembleData}
#'@export
#'@examples
#' N_species <- 4
#' priors <- define_priors(ind_st_var_params = list(25, 0.25),
#'                        ind_st_cor_form = "lkj", #Using an LKJ distribution for individual short-term discrepancies
#'                        ind_st_cor_params = 30, #The parameter is 30
#'                      ind_lt_var_params = list(rep(25,N_species),rep(0.25,N_species)),
#'                      ind_lt_cor_form = "beta",
#'                      ind_lt_cor_params = list(matrix(40,N_species, N_species), matrix(40, N_species, N_species)),
#'                      sha_st_var_exp = 3,
#'                      sha_st_cor_form = "lkj",
#'                      sha_st_cor_params = 30,
#'                    sha_lt_sd = rep(4,N_species))
#'ensemble_data <- EnsembleData(observations = list(SSB_obs, Sigma_obs),
#'                              simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "mizer")),
#'                               priors = priors)
EnsembleData <- function(observations, simulators, priors) {
  # Checks that the input data can be converted into the form required for Stan.
  validate_data(observations, simulators, priors)


  stan_data_sims <- generate_simulator_stan_data(observations, simulators)
  stan_data_priors <- priors[- which(names(priors) == "prior_correlations" )]
  stan_data_priors <- append(stan_data_priors,
                             generate_correlation_priors_stan_data(priors$prior_correlations))

  stan_input <- append(stan_data_sims, stan_data_priors)


  ensemble_data <- new('EnsembleData', stan_input = stan_input,
                       observations = observations,
                       simulators = simulators,
                       priors = priors)
  return(ensemble_data)
}

#### Class definition ####
#' A class to hold the Ensemble data
#'
#' A class that holds the observation data, simulator outputs, and prior information
#' ready to fit the ensemble model in Stan.
#'
#' A new `EnsembleData` object can be created with the [EnsembleData()]
#' constructor, but it is not necessary to do this if using the `fit_ensemble_model`
#' function, as it will be created automatically.
#'
#' There is only one slot for this object, containing a list of data in the required
#' form to fit the Stan model.
#'
#' @slot stan_input A list of parameters in the correct form to fit the ensemble
#' model in Stan.
#' @slot observations A `list` of observations and associated covariances.
#' See the `EnsembleData` constructor details.
#' @slot simulators A `list` of simulator output and covariance pairs. See the `EnsembleData` constructor details.
#' @slot priors A `list` of data specifying the prior distributions for the ensemble.
#' See the `EnsembleData` constructor details.
#'
#' @export
setClass(
  "EnsembleData",
  slots = c(
    stan_input = "list",
    observations = "list",
    simulators = "list",
    priors = "list"
  )
)
