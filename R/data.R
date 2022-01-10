#' Ecopath with EcoSim SSB
#'
#' A dataset containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod,
#' and Sole in the North Sea between 1991-2050 under an MSY fishing scenario from EcoPath with EcoSim.
#'
#' @format A `data frame` with 60 rows and 4 variables:
#' \describe{
#'   \item{`N.Pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references C Walters, V Christensen and D Pauly (1997) Structuring dynamic
#' models of exploited ecosystems from trophic mass-balance assessments
#' Rev. Fish Biol. Fish., 7, pp. 139-172
"SSB_ewe"

#' FishSUMS SSB
#'
#' A `data frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, and Cod
#' in the North Sea between 1984-2050 under an MSY fishing scenario from FishSUMS. Note that FishSUMS does not
#' produce outputs for Sole.
#'
#' @format A data frame with 67 rows and 3 variables:
#' \describe{
#'   \item{`N.Pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#' }
#' @references Speirs, D.C., Guirey, E.J., Gurney, W.S.C. and Heath, M.R. (2010).
#' A length structured partial ecosystem model for cod in the North Sea.
#' Fisheries Research 106, 474-494.
"SSB_fs"

#' LeMans SSB
#'
#' A `data frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod, and Sole
#' in the North Sea between 1986-2050 under an MSY fishing scenario from the LeMaRns package.
#'
#' @format A data frame with 65 rows and 4 variables:
#' \describe{
#'   \item{`N.Pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references   Michael A. Spence, Hayley J. Bannister, Johnathan E. Ball, Paul J. Dolder
#' and Robert B. Thorpe (2019). LeMaRns: Length-Based Multispecies Analysis
#' by Numerical Simulation. R package version 0.1.2.
#' \href{https://CRAN.R-project.org/package=LeMaRns}{CRAN}
"SSB_lm"


#' mizer SSB
#'
#' A `data.frame` containing the predictions for spawning stock biomass of Norway Pout, Herring, Cod, and Sole
#' in the North Sea between 1984-2050 under an MSY fishing scenario from mizer.
#'
#' @format A data frame with 67 rows and 4 variables:
#' \describe{
#'   \item{`N.Pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
#'   \item{`Herring`}{Spawning stock biomass of Herring, in log tonnes.}
#'   \item{`Cod`}{Spawning stock biomass of Cod, in log tonnes.}
#'   \item{`Sole`}{Spawning stock biomass of Sole, in log tonnes.}
#' }
#' @references F. Scott, J.L. Blanchard and K.H. Andersen. mizer: an R package for
#' multispecies, trait-based and community size spectrum ecological
#' modelling. Methods in Ecology and Evolution 5(10) 1121-1125 (2014).
"SSB_miz"

#' Stock assessment SSB
#'
#' A `data frame` containing the single species stock assessment estimates of spawning stock
#' biomass of Norway Pout, Herring, Cod, and Sole in the North Sea between 1984-2017.
#'
#' @format A data frame with 34 rows and 4 variables:
#' \describe{
#'   \item{`N.Pout`}{Spawning stock biomass of Norway Pout, in log tonnes.}
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
"Sigma_ewe"

#' FishSUMS Sigma
#'
#' A 3x3 covariance matrix quantifying the parameter uncertainty of FishSUMS
"Sigma_fs"

#' LeMans Sigma
#'
#' A 4x4 covariance matrix quantifying the parameter uncertainty of LeMans
"Sigma_lm"

#' mizer Sigma
#'
#' A 4x4 covariance matrix quantifying the parameter uncertainty of mizer
"Sigma_miz"

#' Stock assessment Sigma
#'
#' A 4x4 covariance matrix quantifying the covariances of the stock assessment estimates of biomass.
"Sigma_obs"
