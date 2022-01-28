#' Plot the ensemble output
#'
#'Plots the latent truth predicted by the ensemble model, along with simulator outputs and observations.
#'@param sample An `EnsembleSample` object output from `generate_sample`.
#'@param variable The name of the variable / species to plot. This can either be a `character` string in the same
#' form as the observation variable, or an index for the column in the observations data frame.
#'@param ggplot A `logical` which plots the graphs using ggplot if `TRUE` and base R if `FALSE`.
#'@param quantiles A `numeric` vector of length 2 giving the quantiles for which to plot ribbons if doing a full sampling of the ensemble model.
#'@return The `ggplot` object.
#'@export
#'@examples
#'fit_sample <- fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                          simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                          priors = priors)
#'sample <- generate_sample(fit_sample, num_samples = 2000)
#'plot(sample)
#'plot(sample, variable = "Cod", quantiles=c(0.2, 0.8))
#'
#'
#'fit_point <-fit_ensemble_model(observations = list(SSB_obs, Sigma_obs),
#'                               simulators = list(list(SSB_ewe, Sigma_ewe, "EwE"),
#'                                            list(SSB_fs,  Sigma_fs, "FishSUMS"),
#'                                            list(SSB_lm,  Sigma_lm, "LeMans"),
#'                                            list(SSB_miz, Sigma_miz, "Mizer")),
#'                                priors = priors,
#'                                full_sample = FALSE)
#'sample1 <- generate_sample(fit_point, num_samples = 2000)
#'plot(sample1, variable="Herring")
plot.EnsembleSample <- function(sample, variable=1, quantiles=c(0.05, 0.95)){

  fit <- sample@ensemble_fit
  ensemble_data <- fit@ensemble_data
  observations <- ensemble_data@observations[[1]]
  simulators <- ensemble_data@simulators
  stan_input <- ensemble_data@stan_input

  # Find which variable we are interested in
  if(is.double(variable)){
    variable <- colnames(observations)[variable]
  }
  if (is.na(variable)
      || !(variable %in% colnames(observations) )){
    stop(paste0("Invalid variable. This should be the name of a species or an index less than ",
                ncol(observations)))
  }

  #We have to use the stan_input version because `observations` does not have data in the future
  df <- fit@ensemble_data@stan_input$observations[, which(colnames(observations) == variable)]

  #Ensemble outputs
  if(!is.null(fit@samples)){
    df <- cbind(df, apply(sample@mle[, which(colnames(observations) == variable), ], 1, median))
    df <- cbind(df, apply(sample@sample[, which(colnames(observations) == variable), ], 1, quantile, min(quantiles)))
    df <- cbind(df, apply(sample@sample[, which(colnames(observations) == variable), ], 1, quantile, max(quantiles)))
    colnames(df)[(ncol(df) - 2):ncol(df)] <- c("Ensemble true value", "EnsembleLower", "EnsembleUpper")

  }else{
    df <- cbind(df, sample@mle[, which(colnames(observations) == variable)])
    colnames(df)[ncol(df)] <- "Ensemble true value"
  }

  #Simulator outputs
  for (i in 1:length(simulators)) {
    simulator <- simulators[[i]]

    #Skip simulators that dont have the species
    if (!(variable %in% colnames(simulator[[1]]))){
      next
    }

    #Find the right output
    relative_index <- which(colnames(simulator[[1]]) == variable)
    index <- cumsum(stan_input$model_num_species)[i] - stan_input$model_num_species[1] + relative_index
    df <- cbind(df, stan_input$model_outputs[, index])

    #Use the name if available
    colnames(df)[ncol(df)] <- paste0("Simulator ", i)
    if (length(simulator) == 3){
      colnames(df)[ncol(df)] <- simulator[[3]]
    }
  }
  colnames(df)[1] <- "Observations"

  #Years
  start_year <- min(as.double(rownames(observations)))
  df <- cbind(start_year:(start_year+stan_input$time-1), df)
  colnames(df)[1] <- "Year"

  #Keeping zeros makes the graphs nosedive simply because we're missing data/outputs
  #TODO: Worry about legitimate zeros here!
  df <- data.frame(df)
  for (i in 1:ncol(df)) {
    df[df[, i] == 0, i] = NA
  }


  if(!is.null(fit@samples)){
    df <-  reshape2::melt(df, id.vars=c("Year", "EnsembleLower", "EnsembleUpper"), variable.name="Simulator")
    return(plot_values_sample_gg(df, variable))
  }
  df <-  reshape2::melt(df, id.vars=c("Year"), variable.name="Simulator")
  return(plot_values_optimised_gg(df, variable))

}

#Is this needed?
#setMethod("plot", signature(x = "EnsembleSample", y = "missing"),
#          function(x, ...) {
#            plot_values(x, ...)
#          }
#)

plot_values_optimised_gg<- function(df, title){
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x=`Year`, y=`value`, na.rm = TRUE)) +
    ggplot2::geom_line(ggplot2::aes(group=`Simulator`,colour=`Simulator`)) +
    ggplot2::ggtitle(title)
  return(p)
}

plot_values_sample_gg<- function(df, title){
  p <- plot_values_optimised_gg(df, title)

  g <- ggplot2::ggplot_build(p)
  ens_colour <- subset(g$data[[1]], group == 2)
  ens_colour <- unique(ens_colour["colour"])$colour

  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=`EnsembleLower`, ymax =`EnsembleUpper`),
                                alpha=0.2, fill=ens_colour)
  return(p)
}



