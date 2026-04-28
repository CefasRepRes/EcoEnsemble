#' Plot the ensemble output
#'
#'Plots the latent variables predicted by the ensemble model, along with simulator outputs and observations.
#'@param x An `EnsembleSample` object.
#'@param variable The name of the variable to plot. This can either be a `character` string in the same form as the observation variable, or an index for the column in the observations data frame.
#'@param quantiles A `numeric` vector of length 2 giving the quantiles for which to plot ribbons if doing a full sampling of the ensemble model. The default is `c(0.05,0.95)`.
#'@param ... Other arguments passed on to methods. Not currently used.
#'@return The `ggplot` object.
#'@importFrom stats median
#'@importFrom stats quantile
#'@importFrom stats rnorm
#'@importFrom cowplot plot_grid
#'@export
#'@examples
#'\donttest{
#' priors <- EnsemblePrior(4)
#' prior_density <- prior_ensemble_model(priors, M = 4)
#' samples <- sample_prior(observations = list(SSB_obs, Sigma_obs),
#'              simulators = list(list(SSB_miz, Sigma_miz),
#'                                list(SSB_ewe, Sigma_ewe),
#'                                list(SSB_fs, Sigma_fs),
#'                                list(SSB_lm, Sigma_lm)),
#'              priors = priors,
#'              sam_priors = prior_density)
#' plot(samples) #Plot the prior predictive density.
#' plot(samples, variable="Herring")
#'}
plot.EnsembleSample <- function(x, variable = NULL, quantiles=c(0.05, 0.95), ...){
  stan_input <- x@ensemble_fit@ensemble_data@stan_input
  if ("MM" %in% names(stan_input)){
    if(!is.null(variable)){
      return(plot_single_dri(x, variable, quantiles, ...))
    }

    d <- x@ensemble_fit@ensemble_data@priors@d
    plots_all <- lapply(1:d, function(i){plot_single_dri(x, i, quantiles, ...) + ggplot2::theme(legend.position = "none")})
    legend <- cowplot::get_plot_component(plot_single_dri(x, 1), "guide-box-right")
    plots_all <- append(plots_all, list(legend))
  }else{
    if(!is.null(variable)){
      return(plot_single(x, variable, quantiles, ...))
    }

    d <- x@ensemble_fit@ensemble_data@priors@d
    plots_all <- lapply(1:d, function(i){plot_single(x, i, quantiles, ...) + ggplot2::theme(legend.position = "none")})
    legend <- cowplot::get_plot_component(plot_single(x, 1), "guide-box-right")
    plots_all <- append(plots_all, list(legend))
  }
  return(
    do.call(cowplot::plot_grid, plots_all)
  )
}

get_variable <- function(variable, observations){
  ret <- variable
  if(is.double(variable) || is.integer(variable)){
    ret <- colnames(observations)[variable]
    if (abs(variable - round(variable)) > .Machine$double.eps^0.5){
      warning("Non-integer variable specified. variable will be taken as floor(variable). This is done by R, strange isn't it?!")
    }
  }

  if (is.na(ret)
      || !(ret %in% colnames(observations) )){
    stop(paste0("Invalid variable. This should be the name of a variable or an index less than ",
                ncol(observations) + 1))
  }

  return(ret)

}

construct_plot_dataframe <- function(samples, variable, quantiles){

  fit <- samples@ensemble_fit
  ensemble_data <- fit@ensemble_data
  observations <- ensemble_data@observations[[1]]
  simulators <- ensemble_data@simulators
  stan_input <- ensemble_data@stan_input

  #Observations
  df <- tibble::rownames_to_column(observations)[, c("rowname", variable)]
  colnames(df) <- c("Year", "Observations")

  #Simulators
  for (i in 1:length(simulators)) {
    simulator <- simulators[[i]]

    #Skip simulators that dont have the VoI
    if (!(variable %in% colnames(simulator[[1]]))){
      next
    }

    df_sim <- tibble::rownames_to_column(simulator[[1]], var = "Year")[, c("Year", variable)]
    #Use the name if available
    colnames(df_sim)[2] <- paste0("Simulator ", i)
    if (length(simulator) == 3){
      colnames(df_sim)[2] <- simulator[[3]]
    }


    df <- dplyr::full_join(df, df_sim, by = "Year")

  }

  #Ensemble
  var_index = which(colnames(observations) == variable)
  df$Year <- as.numeric(df$Year)
  times <- sort(unique(df$Year))
  if(!is.null(samples@samples)){

    # The ensemble outputs need to be the first columns to ensure they are coloured consistently when some variables are missing.
    df_ensemble <-cbind(times, apply(samples@samples[, var_index, ], 1, median, na.rm = TRUE),
                        apply(samples@samples[, var_index, ], 1, quantile, min(quantiles), na.rm = TRUE),
                        apply(samples@samples[, var_index, ], 1, quantile, max(quantiles), na.rm = TRUE))%>%
                        data.frame()


    colnames(df_ensemble) <- c("Year","Ensemble Model Prediction", "Lower", "Upper")
    df <- merge(df_ensemble,df)

    df <-  reshape2::melt(df, id.vars=c("Year", "Lower", "Upper"), variable.name="Simulator")

    #We have zero-width ribbons for the models and observations to avoid clutterling the scene
    df[df$Simulator != "Ensemble Model Prediction", c("Lower", "Upper")] <- df[df$Simulator != "Ensemble Model Prediction", "value"]

  }else{
    df_ensemble <-data.frame(times, samples@mle[, var_index])

    colnames(df_ensemble) <- c("Year","Ensemble Model Prediction")
    df <- merge(df_ensemble,df)

    df <-  reshape2::melt(df, id.vars=c("Year"), variable.name="Simulator")

  }
  return(df)
}

# Fudge to get past "no visible binding for global variable" in R-CMD check
utils::globalVariables(c("Year", "Lower", "Upper", "value", "Simulator"))

plot_single <- function(samples, variable=1, quantiles=c(0.05, 0.95), ...){

  fit <- samples@ensemble_fit
  ensemble_data <- fit@ensemble_data
  observations <- ensemble_data@observations[[1]]
  simulators <- ensemble_data@simulators
  stan_input <- ensemble_data@stan_input

  variable <- get_variable(variable, observations)

  df <- construct_plot_dataframe(samples, variable, quantiles)

  var_index = which(colnames(observations) == variable)
  if(!is.null(samples@samples)){
    return(plot_values_sample_gg(df, variable))
  }
  return(plot_values_optimised_gg(df, variable))
}

plot_values_optimised_gg<- function(df, title, ...){
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=`Year`, y=`value`, na.rm = TRUE), ...) + ggplot2::geom_line(ggplot2::aes(group=`Simulator`,colour=`Simulator`), na.rm = TRUE) + ggplot2::ggtitle(title)
  return(p)
}

plot_values_sample_gg<- function(df, title, ...){
    return(plot_values_optimised_gg(df, title, ...) +
             ggplot2::geom_ribbon(ggplot2::aes(ymin=`Lower`, ymax =`Upper`, fill = `Simulator`), alpha=0.2,na.rm = TRUE)) +
    ggplot2::theme(legend.position = "right")
}
