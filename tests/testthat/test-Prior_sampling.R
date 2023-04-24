test_that("prior test",{
	testthat::skip_on_cran()
  priors1 <- EnsemblePrior(4,ind_st_params = IndSTPrior(parametrisation_form = "lkj", var_params= list(1,1), cor_params = 10, AR_params = c(2, 2)))
  set.seed(1234)
  prior_density <- prior_ensemble_model(priors1, M = 4)
  samples <- sample_prior(observations = list(SSB_obs, Sigma_obs),
                          simulators = list(list(SSB_miz, Sigma_miz),
                                            list(SSB_ewe, Sigma_ewe),
                                            list(SSB_fs, Sigma_fs),
                                            list(SSB_lm, Sigma_lm)),
                          priors = priors1,
                          sam_priors = prior_density)
  plot(samples)
  priors2 <- EnsemblePrior(4)
  prior_density1 <- prior_ensemble_model(priors1, M = 4,full_sample = FALSE)
  prior_density2 <- prior_ensemble_model(priors2, M = 4,full_sample = TRUE)
  expect_error(prior_ensemble_model(priors2, M = 4,full_sample = FALSE),"It is possible to generate a point estimate for the prior if the individual short-term discrepancy prior follows a hierarchical parameterisation. Please generate a full sample using 'full_sample=TRUE'.")
})

test_that("Ensemble prior test",{
  d <-4
  prior1 <- EnsemblePrior(d=d)
  iLT <- IndLTPrior("lkj",list(1,1),1)
  iST <- IndSTPrior(parametrisation_form = "hierarchical", var_params= list(-3, 1, 8, 4), cor_params = list(0.1, 0.1, 0.1, 0.1), AR_params = c(2, 2))
  sST <- ShaSTPrior(parametrisation_form = "lkj", var_params = list(1, 10), cor_params = 1, AR_params = c(2, 2))
  tp <- TruthPrior(d, initial_mean = 0, initial_var = 100, rw_covariance = list(2* d, diag(d)))
  prior2 <- EnsemblePrior(d=d,ind_st_params = iST,ind_lt_params = iLT,sha_st_params = sST,truth_params = tp)
  expect_equal(prior1,prior2)
})
