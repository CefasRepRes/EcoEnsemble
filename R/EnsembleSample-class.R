#' @rdname EnsembleSample
#' @export
EnsembleSample <- function(ensemble_fit, mle, sample) {
  ensemble_sample <- new('EnsembleSample',
                         ensemble_fit = ensemble_fit,
                         mle = mle,
                         sample = sample)
  return(ensemble_sample)
}

#### Class definition ####
#' A class to hold the samples of the ensemble model from the Kalman filter
#'
#'
#' `EnsembleSample` objects are generated using the `generate_sample` function.
#'
#'
#' @slot ensemble_fit An `EnsembleFit` object containing the fitted ensemble model.
#' @slot mle An `array` of dimension \eqn{T \times (M + 2)\times N_{sample}} containing MLe point estimates from the `ensemble_fit` object, where \eqn{T} is the total time, \eqn{M} is the number of simulators and \eqn{N_{sample}} is the number of samples. For each time step, the `t`th element of the array is a `matrix` where each column is a sample and the rows are the variables:
#' \deqn{\left( y^{(t)}, \eta^{(t)}, z_1^{(t)}, z_2^{(t)}, \ldots, z_M^{(t)}\right)'}
#' where \eqn{y^{(t)}} is the ensemble model's prediction of the latent truth value at time \eqn{t},
#' \eqn{\eta^{(t)}} is the shared short-term discrepancy at time \eqn{t},
#' \eqn{z_i^{(t)}} is the individual short-term discrepancy of simulator \eqn{i} at time \eqn{t}.

#' @slot sample An `array` of dimension \eqn{T \times (M + 2)\times N_{sample}} containing samples from the `ensemble_fit` object, where \eqn{T} is the total time, \eqn{M} is the number of simulators and \eqn{N_{sample}} is the number of samples. For each time step, the `t`th element of the array is a `matrix` where each column is a sample and the rows are the variables:
#' \deqn{\left( y^{(t)}, \eta^{(t)}, z_1^{(t)}, z_2^{(t)}, \ldots, z_M^{(t)}\right)'}
#' where \eqn{y^{(t)}} is the ensemble model's prediction of the latent truth value at time \eqn{t},
#' \eqn{\eta^{(t)}} is the shared short-term discrepancy at time \eqn{t},
#' \eqn{z_i^{(t)}} is the individual short-term discrepancy of simulator \eqn{i} at time \eqn{t}.
#' @seealso \code{\link{EnsembleSample}}, \code{\link{generate_sample}}
#' @rdname EnsembleSample
#' @export
setClass(
  "EnsembleSample",
  slots = c(
    ensemble_fit = "EnsembleFit",
    mle = "array",
    sample = "array"
  )
)
