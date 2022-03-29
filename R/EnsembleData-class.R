#' Constructor for the `EnsembleData` class
#'
#' A constructor for the `EnsembleData` class. This is used to convert input data into the required form for `fit_ensemble_model`.
#'@param observations A `list` of length 2 containing observations and a covariance matrix. The first element is a `data.frame` or `matrix` with each column giving observations of each output of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is is a \eqn{d \times d} `matrix` where \eqn{d} is the number of columns of the observations data frame / matrix. This matrix is the covariance matrix of the observations.
#'@param simulators A `list` with length equal to the number of simulators. For each simulator, there is a `list` of 2 objects containing the simulator output and covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a simulator outputs of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is a \eqn{n_k \times n_k} `matrix` where \eqn{n_k} is the number of columns of the simulators output data frame / matrix. This matrix is the covariance matrix of the simulator outputs.
#'@param  priors An `EnsemblePrior` object specifying the prior distributions for the ensemble.
#'@return An object of class \code{EnsembleData}
#'@seealso \code{\linkS4class{EnsembleData}}, \code{\link{EnsemblePrior}}, \code{\link{fit_ensemble_model}}
#'@export
#'@include EnsemblePrior-class.R
#'@examples
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
  stan_input <- append(stan_data_sims, priors@priors_stan_input)


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
#' @slot stan_input A `list` of parameters in the correct form to fit the ensemble model in Stan.
#' @slot observations A `list` of length 2 containing observations and a covariance matrix. The first element is a `data.frame` or `matrix` with each column giving observations of each output of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the observations.
#' @slot simulators A `list` with length equal to the number of simulators. For each simulator, there is a `list` of 2 objects containing the simulator output and covariance matrix. The first element is a `data.frame` or `matrix` with each column giving a simulator outputs of interest and each row a time. Rows should be named with the times and columns should be named the variables. The second element is the covariance matrix of the simulator outputs.
#' @slot priors An `EnsemblePrior` object specifying the prior distributions for the ensemble.
#' @seealso \code{\link{EnsembleData}}, \code{\link{EnsemblePrior}}, \code{\link{fit_ensemble_model}}
#' @export
setClass(
  "EnsembleData",
  slots = c(
    stan_input = "list",
    observations = "list",
    simulators = "list",
    priors = "EnsemblePrior"
  )
)
