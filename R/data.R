#' Ecopath with EcoSim SSB
#'
#' A dataset containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod, and Sole in the North Sea between 1991-2050 under an MSY fishing scenario from EcoPath with EcoSim.
#'
#' @format A `data.frame` with 60 rows and 4 variables:
#' \describe{
#'   \item{`N.pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references ICES (2016). Working Group on Multispecies Assessment Methods (WGSAM). Technical report, International Council for Exploration of the Seas.
"SSB_ewe"

#' FishSUMS SSB
#'
#' A `data frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, and Cod
#' in the North Sea between 1984-2050 under an MSY fishing scenario from FishSUMS. Note that FishSUMS does not
#' produce outputs for Sole.
#'
#' @format A `data.frame` with 67 rows and 3 variables:
#' \describe{
#'   \item{`N.pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#' }
#' @references Speirs, D., Greenstreet, S., and Heath, M. (2016). Modelling the effects of fishing on the North Sea fish community size composition. Ecological Modelling, 321, 35–45
"SSB_fs"

#' LeMans SSB
#'
#' A `data.frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod, and Sole
#' in the North Sea between 1986-2050 under an MSY fishing scenario from the LeMaRns package.
#'
#' @format A `data.frame` with 65 rows and 4 variables:
#' \describe{
#'   \item{`N.pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references Thorpe, R. B., Le Quesne, W. J. F., Luxford, F., Collie, J. S., and Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. Methods in Ecology and Evolution, 6(1), 49–58.
"SSB_lm"


#' mizer SSB
#'
#' A `data.frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod, and Sole
#' in the North Sea between 1984-2050 under an MSY fishing scenario from mizer.
#'
#' @format A data frame with 67 rows and 4 variables:
#' \describe{
#'   \item{`N.pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references Blanchard, J. L., Andersen, K. H., Scott, F., Hintzen, N. T., Piet, G., and Jennings, S. (2014). Evaluating targets and trade-offs among fisheries and conservation objectives using a multispecies size spectrum model. Journal of Applied Ecology, 51(3), 612–622
"SSB_miz"

#' Stock assessment SSB
#'
#' A `data.frame` containing the single species stock assessment estimates of spawning stock
#' biomass of Norway Pout, Herring, Cod, and Sole in the North Sea between 1984-2017.
#'
#' @format A data frame with 34 rows and 4 variables:
#' \describe{
#'   \item{`N.pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references Herring Assessment Working Group for the Area South of 62 N (HAWG).Technical report,
#' ICES Scientific Reports. ACOM:07. 960 pp, ICES, Copenhagen.
#' @references Report of the Working Group on the Assessment of Demersal Stocks inthe North Sea and
#' Skagerrak. Technical report, ICES Scientific Reports. ACOM:22.pp, ICES, Copenhagen.
"SSB_obs"


#' Ecopath with EcoSim Sigma
#'
#' A 4x4 covariance matrix quantifying the parameter uncertainty of Ecopath with EcoSim
#' @format A 4x4 `matrix`.
#' @references Mackinson, S., Platts, M., Garcia, C., and Lynam, C. (2018). Evaluating the fishery and ecological consequences of the proposed North Sea multi-annual plan. PLOS ONE, 13(1), 1–23.
"Sigma_ewe"

#' FishSUMS Sigma
#'
#' A 3x3 covariance matrix quantifying the parameter uncertainty of FishSUMS
#' @format A 3x3 `matrix`.
#' @references Spence, M. A., Blanchard, J. L., Rossberg, A. G., Heath, M. R., Heymans, J. J., Mackinson, S., Serpetti, N., Speirs, D. C., Thorpe, R. B., and Blackwell, P. G. (2018). A general framework for combining ecosystem models. Fish and Fisheries, 19(6), 1031–1042.
"Sigma_fs"

#' LeMans Sigma
#'
#' @format A 4x4 `matrix`.
#' A 4x4 covariance matrix quantifying the parameter uncertainty of LeMans
#' @references   Thorpe, R. B., Le Quesne, W. J. F., Luxford, F., Collie, J. S., and Jennings, S. (2015). Evaluation and management implications of uncertainty in a multispecies size-structured model of population and community responses to fishing. Methods in Ecology and Evolution, 6(1), 49–58.
"Sigma_lm"

#' mizer Sigma
#'
#' @format A 4x4 `matrix`.
#' A 4x4 covariance matrix quantifying the parameter uncertainty of mizer
#' @references Spence, M. A., Blackwell, P. G., and Blanchard, J. L. (2016). Parameter uncertainty of a dynamic multispecies size spectrum model. Canadian Journal of Fisheries and Aquatic Sciences, 73(4), 589–59

"Sigma_miz"

#' Stock assessment Sigma
#'
#' A 4x4 covariance matrix quantifying the covariances of the stock assessment estimates of biomass.
#' @format A 4x4 `matrix`.
#' @references Herring Assessment Working Group for the Area South of 62 N (HAWG).Technical report,
#' ICES Scientific Reports. ACOM:07. 960 pp, ICES, Copenhagen.
#' @references Report of the Working Group on the Assessment of Demersal Stocks inthe North Sea and
#' Skagerrak. Technical report, ICES Scientific Reports. ACOM:22.pp, ICES, Copenhagen.

"Sigma_obs"
