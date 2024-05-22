plot_withdrivers.EnsembleSample <- function(x, variable = NULL, quantiles=c(0.05, 0.95), drivers = TRUE, ...){
  if (drivers == FALSE) {
    if(!is.null(variable)){
      return(plot_single(x, variable, quantiles, ...))
    }

    d <- x@ensemble_fit@ensemble_data@priors@d
    plots_all <- lapply(1:d, function(i){plot_single(x, i, quantiles, ...) + ggplot2::theme(legend.position = "none")})
    legend <- cowplot::get_plot_component(plot_single(x, 1), "guide-box-right")
    plots_all <- append(plots_all, list(legend))
  }else{
    if(!is.null(variable)){
      return(plot_single_dri(x, variable, quantiles, ...))
    }

    d <- x@ensemble_fit@ensemble_data@priors@d
    plots_all <- lapply(1:d, function(i){plot_single_dri(x, i, quantiles, ...) + ggplot2::theme(legend.position = "none")})
    legend <- cowplot::get_plot_component(plot_single(x, 1), "guide-box-right")
    plots_all <- append(plots_all, list(legend))
  }
  return(
    do.call(cowplot::plot_grid, plots_all)
  )
}


construct_plot_dataframe_dri <- function(samples, variable, quantiles){

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
      for (j in 1:length(simulator[[1]])){
        driver <- simulator[[1]][[j]]
        #Skip simulators that don't have the VoI
        if (!(variable %in% colnames(driver))){
          next
        }

        df_sim <- tibble::rownames_to_column(driver, var = "Year")[, c("Year", variable)]
        #Use the name if available
        colnames(df_sim)[2] <- paste0("Simulator ", i, "+", "Driver ", j)
        if (length(simulator) >= 3){
          if (length(simulator) == 4){
            colnames(df_sim)[2] <- paste0(simulator[[3]], "+", simulator[[4]][[j]])
          } else {
            colnames(df_sim)[2] <- paste0(simulator[[3]], "+", "Driver ", j)
          }
        }
        df <- dplyr::full_join(df, df_sim, by = "Year")
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

plot_single_dri <- function(samples, variable=1, quantiles=c(0.05, 0.95), ...){

  fit <- samples@ensemble_fit
  ensemble_data <- fit@ensemble_data
  observations <- ensemble_data@observations[[1]]
  simulators <- ensemble_data@simulators
  stan_input <- ensemble_data@stan_input

  variable <- get_variable(variable, observations)

  df <- construct_plot_dataframe_dri(samples, variable, quantiles)

  var_index = which(colnames(observations) == variable)
  if(!is.null(samples@samples)){
    return(plot_values_sample_gg(df, variable))
  }
  return(plot_values_optimised_gg(df, variable))
}


