test_that("get ensemble model works", {

  priors_hier <- EnsemblePrior(4)
  priors_lkj <- EnsemblePrior(4, ind_st_params = IndSTPrior("lkj", list(3, 2), 3, AR_params = c(2,4)))

  expect_equal(get_mcmc_ensemble_model(priors_hier)@model_name, "ensemble_model_hierarchical")
  expect_equal(get_mcmc_ensemble_model(priors_hier, likelihood = FALSE)@model_name, "ensemble_prior_hierarchical")
  expect_equal(get_mcmc_ensemble_model(priors_lkj)@model_name, "ensemble_model")
  expect_equal(get_mcmc_ensemble_model(priors_lkj, likelihood = FALSE)@model_name, "ensemble_prior")

  expect_equal(get_mcmc_ensemble_model(priors_hier, drivers = TRUE)@model_name, "ensemble_model_hierarchical_withdrivers")
  expect_equal(get_mcmc_ensemble_model(priors_hier, likelihood = FALSE, drivers = TRUE)@model_name, "ensemble_prior_hierarchical_withdrivers")
  expect_equal(get_mcmc_ensemble_model(priors_lkj, drivers = TRUE)@model_name, "ensemble_model_withdrivers")
  expect_equal(get_mcmc_ensemble_model(priors_lkj, likelihood = FALSE, drivers = TRUE)@model_name, "ensemble_prior_withdrivers")
})
