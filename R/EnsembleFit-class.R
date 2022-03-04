#' @rdname EnsembleFit
#' @export
EnsembleFit <- function(ensemble_data, samples = NULL, point_estimate = NULL) {
  ensemble_data <- new('EnsembleFit',
                       ensemble_data = ensemble_data,
                       samples = samples,
                       point_estimate = point_estimate)
  return(ensemble_data)
}


setClassUnion("stanfit_or_null", c("stanfit", "NULL"))
setClassUnion("list_or_null", c("list", "NULL"))


#### Class definition ####
#' A class to hold the samples or point estimates from the ensemble model.
#'
#' An `EnsembleFit` object is returned by the `fit_ensemble_model` function.
#' The object contains a slot for the `EnsembleData` object originally used to
#' fit the ensemble model. The `samples` slot contains the samples from the MCMC
#' if a full sampling was completed, otherwise the `point_estimate` slot contains
#' information about a point estimate.
#'
#'
#' @slot ensemble_data An `EnsembleData` object encapsulating the data used to fit the ensemble model.
#' @slot samples A `stanfit` object containing the samples drawn from the fitted model. The default value is `NULL`.
#' @slot point_estimate A `list` output of the optimised model. The default value is `NULL`.
#'
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#' @seealso \code{\link{EnsembleFit}}, \code{\link{fit_ensemble_model}}, \code{\link{generate_sample}}
#' @rdname EnsembleFit
#' @export
setClass(
  "EnsembleFit",
  slots = c(
    ensemble_data = "EnsembleData",
    samples = "stanfit_or_null",
    point_estimate = "list_or_null"
  )
)
