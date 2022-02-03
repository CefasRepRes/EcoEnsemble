#' Constructor for the `EnsembleData` class
#'
#' A constructor for the `EnsembleData` class. This is used to convert input data into the required form to fit the ensemble model.
#'@param observations A `list` 2 objectes containing observations and a covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a observations of each output of interest and each row a different time. Rows should be named with the times and columns should be named the appropriate variables. The second element is the covariance matrix of the observations.
#'@param simulators A `list` of length the number of simulators. For each simulator, there is a `list` of 2 objects containing simulator output and covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a simulator outputs of interest and each row a different time. Rows should be named with the times and columns should be named the appropriate variables. The second element is the covariance matrix of the simulator outputs.
#'@param priors A `list` of data specifying the prior distributions for the ensemble.
#'See the `define_priors` function for the required form.
#'@details Rows should be named with the times and columns should be named the appropriate variables. It is fine to have missing years or outputs, however for any simulator output there should be at least one observation in the observation data frame or matrix. Missing outputs and years are encoded by simply not including the corresponding row or column in the data frame or matrix.
#'@return An object of class \code{\linkS4class{EnsembleData}}
#'@seealso \code{\linkS4class{EnsembleData}}, \code{\link{define_priors}}
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
#' A class that holds the observation data, simulator outputs, and prior information ready to fit the ensemble model in Stan.
#'
#' @slot stan_input A list of parameters in the correct form to fit the ensemble
#' model in Stan.
#' @slot observations A `list` of observations and associated covariances. --- same as above
#' @slot simulators A `list` of simulator output and covariance pairs --- same as above
#' @slot priors A `list` of data specifying the prior distributions for the ensemble --- same as above
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
