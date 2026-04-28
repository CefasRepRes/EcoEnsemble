test_that("prior combinations",{
  d<- 2
  M <- 2
  Time <- 10
  Times_obs <- round(Time * 0.8)
  obs <- matrix(c(0.57,0.75,0.9,1,1.07,1.09,1.04,0.94,0.99,0.54,0.57,1.09,1.67,2.01,2.17,2.33),ncol=d)
  obs.cov <- matrix(c(0.07,-0.02,-0.02,0.01),ncol=d)
  model.cov <- array(0,dim=c(M,d,d))
  models_output <- array(NA,dim=c(M,Time,d))
  models_output[1,,] <- c(-0.14,0.3,0.78,1.19,1.38,1.27,1.04,0.97,1.27,1.85,3.75,3.18,3.01,3.33,3.91,4.41,4.71,4.97,5.42,6.08)
  models_output[2,,] <- c(-0.94,-0.55,-0.25,-0.15,-0.26,-0.49,-0.7,-0.73,-0.52,-0.17,1.48,0.93,0.72,0.98,1.56,2.14,2.54,2.78,2.99,3.25)
  model.cov[1,,] <- c(0.06,0,0,0.03)
  model.cov[2,,] <- c(0.04,0.01,0.01,0.02)

  val_obs <- data.frame(obs); cov_obs <- obs.cov
  val_model_1 <- data.frame(models_output[1,,]); cov_model_1 <- model.cov[1,,]
  val_model_2 <- data.frame(models_output[2,,]); cov_model_2 <- model.cov[2,,]


  #Set the dimnames to ensure EcoEnsemble can identify the information.
  SPECIES_NAMES <- c("Species 1", "Species 2")
  dimnames(val_obs) <- list(paste(1:Times_obs), SPECIES_NAMES)
  dimnames(val_model_1) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_2) <- list(paste(1:Time), SPECIES_NAMES)


  dimnames(cov_obs) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_1) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_2) <- list(SPECIES_NAMES, SPECIES_NAMES)


  priors1 <- EnsemblePrior(2,
                           ind_st_params = IndSTPrior(parametrisation_form = "inv_wishart", var_params= list(1,1), cor_params = list(10,diag(2)), AR_params = c(2, 2)),
                           ind_lt_params = IndLTPrior("beta",list(10,5),list(matrix(5, 2, 2),matrix(2, 2, 2))
                           ),
                           sha_st_params = ShaSTPrior("inv_wishart",list(2, 1/3),list(10, diag(2))))

  set.seed(1234)
  samples <- sample_prior(observations = list(val_obs, cov_obs),
                          simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                            list(val_model_2, cov_model_2, "Model 2")
                          ),
                          priors = priors1,
                          full_sample = FALSE)
  fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                            simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                              list(val_model_2, cov_model_2, "Model 2")
                            ),
                            priors = priors1,
                            full_sample = FALSE)


  priors2 <- EnsemblePrior(2,
                           ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(c(1,2),c(1,1)), cor_params = list(matrix(50, 2, 2),matrix(50, 2, 2))
                                                      , AR_params = c(2, 2)),
                           ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
                           ),
                           sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))
                           ))

  for (sampler in c("kalman", "explicit")){

    samples <- sample_prior(observations = list(val_obs, cov_obs),
                            simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                              list(val_model_2, cov_model_2, "Model 2")
                            ),
                            priors = priors2,
                            full_sample = FALSE)
    fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                              simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                list(val_model_2, cov_model_2, "Model 2")
                              ),
                              priors = priors2, sampler = sampler,
                              full_sample = FALSE)
    priors3 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = 10, AR_params = c(2, 2)))
    fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                              simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                list(val_model_2, cov_model_2, "Model 2")
                              ),
                              priors = priors3, sampler = sampler,
                              full_sample = FALSE)

    #Hierarchical priors
    priors4 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical", var_params= list(7,4,6,5), cor_params = list(0.4, 0.2, 0.7, 0.1), AR_params = c(2, 2)))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                                 list(val_model_2, cov_model_2, "Model 2")
                                               ),
                                               priors = priors4,sampler = sampler,
                                               control = list(adapt_delta = 0.9), full_sample = TRUE, chains = 1, iter = 4))
    priors5 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical_beta_conjugate", var_params= list(7,4,6,5), cor_params = list(0.9, 0.9, 0.5), AR_params = c(2, 2)))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                                 list(val_model_2, cov_model_2, "Model 2")
                                               ),
                                               priors = priors5, sampler = sampler,
                                               control = list(adapt_delta = 0.9), full_sample = TRUE, chains = 1, iter = 4))


    #Errors
    expect_error(EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(1,1), cor_params = list(10,diag(2)), AR_params = c(2, 2))),"Invalid beta parameters for individual short-term correlation matrix priors. These should be square, symmetric matrices with the same dimension as the number of model outputs.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "inv_wishart", var_params= list(1,1), cor_params = list(10), AR_params = c(2, 2))),"Invalid inverse Wishart parameters for individual short-term correlation matrix prior. This should be a numeric giving the degrees of freedom and a d x d scale matrix, with degrees of freedom > d-1.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid lkj parameter for individual short-term correlation matrix prior. This should be a list with exactly one real number.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid lkj parameter for individual short-term correlation matrix prior. This should be a list with exactly one real number.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parameters for individual short-term correlation matrix priors. This should be a list of length 4.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical_beta_conjugate", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parameters for individual short-term correlation matrix priors. This should be a list of length 3.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "boro", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parametrisation choice for priors. Prior parametrisation forms should be one of 'lkj', 'inv_wishart', 'beta', 'hierarchical' or 'hierarchical_beta_conjugate'. The 'hierarchical' and 'hierarchical_beta_conjugate' options are only available for individual short-term discrepancies. Prior choice: boro")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = 10, AR_params = c(2, 2,2))),"Invalid autoregressive parameters specified. This should be a numeric of length 2 giving the shape parameters of the Beta distribution used as a prior for the AR parameters.")
    expect_error (EcoEnsemble:::validate_data(observations = list(val_obs, cov_obs),
                                              simulators = list(val_model_1, cov_model_1,"Model 1"
                                              ),
                                              priors = priors1))
    val_model_1a <- val_model_1
    names(val_model_1a)[2] <- "Car"
    expect_error(EcoEnsemble:::validate_data(observations = list(val_obs, cov_obs),
                                             simulators = list(list(val_model_1a,cov_model_1,"Model 1"
                                             )),
                                             priors = priors1))
    #Sampling for non-hierarchical options
    priors1 <- EnsemblePrior(2,
                             ind_lt_params = IndLTPrior("beta",list(10,5),list(matrix(5, 2, 2),matrix(2, 2, 2))
                             ),
                             sha_st_params = ShaSTPrior("inv_wishart",list(2, 1/3),list(10, diag(2))))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                                 list(val_model_2, cov_model_2, "Model 2")
                                               ),
                                               priors = priors1, sampler = sampler,
                                               control = list(adapt_delta = 0.9),chains=1,iter=4))
    priors1 <- EnsemblePrior(2,
                             ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
                             ),
                             sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Model 1"),
                                                                 list(val_model_2, cov_model_2, "Model 2")
                                               ),
                                               priors = priors1,sampler = sampler,
                                               control = list(adapt_delta = 0.9),chains=1,iter=4))
  }
})



test_that("prior combinations with drivers",{
  d<- 2
  M <- 2
  MM <- 2
  Time <- 10
  Times_obs <- round(Time * 0.8)
  obs <- matrix(c(0.57,0.75,0.9,1,1.07,1.09,1.04,0.94,0.99,0.54,0.57,1.09,1.67,2.01,2.17,2.33),ncol=d)
  obs.cov <- matrix(c(0.07,-0.02,-0.02,0.01),ncol=d)
  model.cov <- array(0,dim=c(M * MM,d,d))
  models_output <- array(NA,dim=c(M * MM,Time,d))
  models_output[1,,] <- c(-0.14,0.3,0.78,1.19,1.38,1.27,1.04,0.97,1.27,1.85,3.75,3.18,3.01,3.33,3.91,4.41,4.71,4.97,5.42,6.08)
  models_output[2,,] <- c(-0.94,-0.55,-0.25,-0.15,-0.26,-0.49,-0.7,-0.73,-0.52,-0.17,1.48,0.93,0.72,0.98,1.56,2.14,2.54,2.78,2.99,3.25)
  models_output[3,,] <- c(-0.05, 0.35, 0.84, 1.15, 1.18, 1.17, 1.08, 1.16, 1.94, 2.39, 3.93, 3.19, 2.88, 3.17, 3.69, 3.89, 3.99, 4.85, 5.6, 6.2)
  models_output[4,,] <- c(-1.16,-0.75,-0.65,-0.56,-0.66,-0.79,-0.9,-0.43,-0.02, 0.17,1.18,0.27,0.12,0.38,0.46,0.54,0.94,2.08,2.71,3.35)
  model.cov[1,,] <- c(0.06,0,0,0.03)
  model.cov[2,,] <- c(0.04,0.01,0.01,0.02)
  model.cov[3,,] <- c(0.617, 0.011, 0.011, 0.343)
  model.cov[4,,] <- c(0.19, 0.072, 0.072, 0.263)

  val_obs <- data.frame(obs); cov_obs <- obs.cov
  val_model_11 <- data.frame(models_output[1,,]); cov_model_11 <- model.cov[1,,]
  val_model_21 <- data.frame(models_output[2,,]); cov_model_21 <- model.cov[2,,]
  val_model_12 <- data.frame(models_output[3,,]); cov_model_12 <- model.cov[3,,]
  val_model_22 <- data.frame(models_output[4,,]); cov_model_22 <- model.cov[4,,]

  #Set the dimnames to ensure EcoEnsemble can identify the information.
  SPECIES_NAMES <- c("Species 1", "Species 2")
  dimnames(val_obs) <- list(paste(1:Times_obs), SPECIES_NAMES)
  dimnames(val_model_11) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_21) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_12) <- list(paste(1:Time), SPECIES_NAMES)
  dimnames(val_model_22) <- list(paste(1:Time), SPECIES_NAMES)

  dimnames(cov_obs) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_11) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_21) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_12) <- list(SPECIES_NAMES, SPECIES_NAMES)
  dimnames(cov_model_22) <- list(SPECIES_NAMES, SPECIES_NAMES)

  val_model_1 <- list(val_model_11, val_model_12)
  val_model_2 <- list(val_model_21, val_model_22)
  cov_model_1 <- list(cov_model_11, cov_model_12)
  cov_model_2 <- list(cov_model_21, cov_model_22)


  for (sampler in c("explicit", "kalman")){

    priors1 <- EnsemblePrior(2,
                             ind_st_params = IndSTPrior(parametrisation_form = "inv_wishart", var_params= list(1,1), cor_params = list(10,diag(2)), AR_params = c(2, 2)),
                             ind_lt_params = IndLTPrior("beta",list(10,5),list(matrix(5, 2, 2),matrix(2, 2, 2))
                             ),
                             sha_st_params = ShaSTPrior("inv_wishart",list(2, 1/3),list(10, diag(2))))

    set.seed(5678)
    samples <- sample_prior(observations = list(val_obs, cov_obs),
                            simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                              list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                            ),
                            priors = priors1,
                            full_sample = FALSE, drivers = TRUE)

    fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                              simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                              ),
                              priors = priors1,sampler = sampler,
                              full_sample = FALSE, drivers = TRUE)


    priors2 <- EnsemblePrior(2,
                             ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(c(1,2),c(1,1)), cor_params = list(matrix(50, 2, 2),matrix(50, 2, 2))
                                                        , AR_params = c(2, 2)),
                             ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
                             ),
                             sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))
                             ))

    samples <- sample_prior(observations = list(val_obs, cov_obs),
                            simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                              list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                            ),
                            priors = priors2,
                            full_sample = FALSE, drivers = TRUE)
    fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                              simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                              ),
                              priors = priors2,sampler = sampler,
                              full_sample = FALSE, drivers = TRUE)
    priors3 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = 10, AR_params = c(2, 2)))
    samples <- sample_prior(observations = list(val_obs, cov_obs),
                            simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                              list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                            ),
                            priors = priors3,
                            full_sample = FALSE, drivers = TRUE)
    fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                              simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                              ),
                              priors = priors3,sampler = sampler,
                              full_sample = FALSE, drivers = TRUE)

    #Hierarchical priors
    priors4 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical", var_params= list(7,4,6,5), cor_params = list(0.4, 0.2, 0.7, 0.1), AR_params = c(2, 2)))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                                 list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                                               ),
                                               priors = priors4,sampler = sampler,
                                               control = list(adapt_delta = 0.9), full_sample = TRUE, drivers = TRUE, chains = 1, iter = 4))
    priors5 <- EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical_beta_conjugate", var_params= list(7,4,6,5), cor_params = list(0.9, 0.9, 0.5), AR_params = c(2, 2)))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                                 list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                                               ),
                                               priors = priors5,sampler = sampler,
                                               control = list(adapt_delta = 0.9), full_sample = TRUE, drivers = TRUE, chains = 1, iter = 4))


    #Errors
    expect_error(EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "beta", var_params= list(1,1), cor_params = list(10,diag(2)), AR_params = c(2, 2))),"Invalid beta parameters for individual short-term correlation matrix priors. These should be square, symmetric matrices with the same dimension as the number of model outputs.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "inv_wishart", var_params= list(1,1), cor_params = list(10), AR_params = c(2, 2))),"Invalid inverse Wishart parameters for individual short-term correlation matrix prior. This should be a numeric giving the degrees of freedom and a d x d scale matrix, with degrees of freedom > d-1.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid lkj parameter for individual short-term correlation matrix prior. This should be a list with exactly one real number.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid lkj parameter for individual short-term correlation matrix prior. This should be a list with exactly one real number.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parameters for individual short-term correlation matrix priors. This should be a list of length 4.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "hierarchical_beta_conjugate", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parameters for individual short-term correlation matrix priors. This should be a list of length 3.")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "boro", var_params= list(1,1), cor_params = list(10,2), AR_params = c(2, 2))),"Invalid parametrisation choice for priors. Prior parametrisation forms should be one of 'lkj', 'inv_wishart', 'beta', 'hierarchical' or 'hierarchical_beta_conjugate'. The 'hierarchical' and 'hierarchical_beta_conjugate' options are only available for individual short-term discrepancies. Prior choice: boro")
    expect_error (EnsemblePrior(2,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = 10, AR_params = c(2, 2,2))),"Invalid autoregressive parameters specified. This should be a numeric of length 2 giving the shape parameters of the Beta distribution used as a prior for the AR parameters.")
    expect_error (EcoEnsemble:::validate_data(observations = list(val_obs, cov_obs),
                                              simulators = list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")
                                              ),
                                              priors = priors1))
    val_model_1a <- val_model_1
    names(val_model_1a[[1]])[2] <- "Car"
    expect_error(EcoEnsemble:::validate_data(observations = list(val_obs, cov_obs),
                                             simulators = list(list(val_model_1a,cov_model_1,"Simulator 1", c("Driver 1", "Driver 2")
                                             )),
                                             priors = priors1))
    #Sampling for non-hierarchical options
    priors1 <- EnsemblePrior(2,
                             ind_lt_params = IndLTPrior("beta",list(10,5),list(matrix(5, 2, 2),matrix(2, 2, 2))
                             ),
                             sha_st_params = ShaSTPrior("inv_wishart",list(2, 1/3),list(10, diag(2))))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                                 list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                                               ),
                                               priors = priors1,sampler = sampler,
                                               control = list(adapt_delta = 0.9),chains=1,iter=4,drivers=TRUE))
    priors1 <- EnsemblePrior(2,
                             ind_lt_params = IndLTPrior("inv_wishart",list(c(1,2),c(1,1)),list(10,diag(2))
                             ),
                             sha_st_params = ShaSTPrior("beta",list(c(1,2),c(1,1)),list(matrix(5, 2, 2),matrix(2, 2, 2))))
    suppressWarnings(fit <- fit_ensemble_model(observations = list(val_obs, cov_obs),
                                               simulators = list(list(val_model_1, cov_model_1, "Simulator 1", c("Driver 1", "Driver 2")),
                                                                 list(val_model_2, cov_model_2, "Simulator 2", c("Driver 1", "Driver 2"))
                                               ),
                                               priors = priors1,sampler = sampler,
                                               control = list(adapt_delta = 0.9),chains=1,iter=4,drivers=TRUE))
  }
})
