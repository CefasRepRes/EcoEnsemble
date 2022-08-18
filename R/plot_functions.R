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
#'fit_sample <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors)
#'samples <- generate_sample(fit_sample, num_samples = 2000)
#'plot(samples)
#'plot(samples, variable = "Cod", quantiles=c(0.2, 0.8))
#'
#'
#'fit_point <-fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                               simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                                priors = priors,
#'                                full_sample = FALSE)
#'samples1 <- generate_sample(fit_point, num_samples = 2000)
#'plot(samples1, variable="Herring")
#'}
plot.EnsembleSample <- function(x, variable = NULL, quantiles=c(0.05, 0.95), ...){
  if(!is.null(variable)){
    return(plot_single(x, variable, quantiles, ...))
  }

  d <- x@ensemble_fit@ensemble_data@priors@d
  plots_all <- lapply(1:d, function(i){plot_single(x, i, quantiles, ...) + ggplot2::theme(legend.position = "none")})
  legend <- cowplot::get_legend(plot_single(x, 1))
  plots_all <- append(plots_all, list(legend))
  return(
    do.call(cowplot::plot_grid, plots_all)
  )
}

# Fudge to get past "no visible binding for global variable" in R-CMD check
utils::globalVariables(c("Year", "EnsembleLower", "EnsembleUpper", "value", "Simulator"))

plot_single <- function(samples, variable=1, quantiles=c(0.05, 0.95), ...){

  fit <- samples@ensemble_fit
  ensemble_data <- fit@ensemble_data
  observations <- ensemble_data@observations[[1]]
  simulators <- ensemble_data@simulators
  stan_input <- ensemble_data@stan_input

  # Find which variable we are interested in
  if(is.double(variable) || is.integer(variable)){
    if (abs(variable - round(variable)) > .Machine$double.eps^0.5){
      warning("Non-integer variable specified. variable will be taken as floor(variable). This is done by R, strange isn't it?!")
    }
    variable <- colnames(observations)[variable]
  }
  if (is.na(variable)
      || !(variable %in% colnames(observations) )){
    stop(paste0("Invalid variable. This should be the name of a variable or an index less than ",
                ncol(observations) + 1))
  }


  df <- tibble::rownames_to_column(observations)[, c("rowname", variable)]
  colnames(df) <- c("Year", "Observations")
  for (i in 1:length(simulators)) {
    simulator <- simulators[[i]]

    #Skip simulators that dont have the species
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
  var_index = which(colnames(observations) == variable)
  if(!is.null(fit@samples)){
    df <- cbind(df, apply(samples@mle[, var_index, ], 1, median))
    df <- cbind(df, apply(samples@samples[, var_index, ], 1, quantile, min(quantiles), na.rm = TRUE))
    df <- cbind(df, apply(samples@samples[, var_index, ], 1, quantile, max(quantiles), na.rm = TRUE))
    df <- data.frame(df)
    df$Year <- as.numeric(df$Year)
    colnames(df)[(ncol(df) - 2):ncol(df)] <- c("Ensemble Model Prediction", "EnsembleLower", "EnsembleUpper")

    df_plot <-  reshape2::melt(df, id.vars=c("Year", "EnsembleLower", "EnsembleUpper"), variable.name="Simulator")

    df_plot[df_plot$Simulator != "Ensemble Model Prediction", c("EnsembleLower", "EnsembleUpper")] <- c(NA, NA)
    return(plot_values_sample_gg(df_plot, variable, ...))

  }else{
    df <- cbind(df, samples@mle[, var_index])
    df <- data.frame(df)
    df$Year <- as.numeric(df$Year)
    colnames(df)[ncol(df)] <- "Ensemble Model Prediction"

    df_plot <-  reshape2::melt(df, id.vars=c("Year"), variable.name="Simulator")
    return(plot_values_optimised_gg(df_plot, variable, ...))
  }


}

plot_values_optimised_gg<- function(df, title, ...){
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=`Year`, y=`value`, na.rm = TRUE), ...) + ggplot2::geom_line(ggplot2::aes(group=`Simulator`,colour=`Simulator`), na.rm = TRUE) + ggplot2::ggtitle(title)
  return(p)
}

plot_values_sample_gg<- function(df, title, ...){
    return(plot_values_optimised_gg(df, title, ...) + ggplot2::geom_ribbon(ggplot2::aes(ymin=`EnsembleLower`, ymax =`EnsembleUpper`, fill = `Simulator`), alpha=0.2))
}
