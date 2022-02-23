#' Constructor for the `EnsembleData` class
#'
#' A constructor for the `EnsembleData` class. This is used to convert input data into the required form for `fit_ensemble_model`.
#'@param observations A `list` of length 2 containing observations and a covariance matrix. The first element is a `data.frame` or `matrix` with each column giving observations of each output of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the observations.
#'@param simulators A `list` with length equal to the number of simulators. For each simulator, there is a `list` of 2 objects containing the simulator output and covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a simulator outputs of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the simulator outputs.
#'@param priors A `list` of data specifying the prior distributions for the ensemble.
#'See the `define_priors` function for the required form.
#'@return An object of class \code{EnsembleData}
#'@seealso \code{\linkS4class{EnsembleData}}, \code{\link{define_priors}}, \code{\link{fit_ensemble_model}}
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
#' A class that holds the observation data, simulator outputs, and prior information to convert into the required form for `fit_ensemble_model`.
#'
#' @slot stan_input A list of parameters in the correct form to fit the ensemble model in Stan.
#' @slot observations A `list` of length 2 containing observations and a covariance matrix. The first element is a `data.frame` or `matrix` with each column giving observations of each output of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the observations.
#' @slot simulators A `list` with length equal to the number of simulators. For each simulator, there is a `list` of 2 objects containing the simulator output and covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a simulator outputs of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the simulator outputs.
#' @slot priors A `list` of data specifying the prior distributions for the ensemble.
#'See the `define_priors` function for the required form.
#' @seealso \code{\link{EnsembleData}}, \code{\link{define_priors}}, \code{\link{fit_ensemble_model}}
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
