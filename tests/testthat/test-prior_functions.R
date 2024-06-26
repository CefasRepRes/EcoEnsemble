test_that("Correlation prior strings are identified", {
  expect_equal(correlation_form_prior("LKJ"), 0);
  expect_equal(correlation_form_prior("inv_wishart"), 1);
  expect_equal(correlation_form_prior("beta"), 2);
  expect_equal(correlation_form_prior("hierarchical"), 3);
  expect_equal(correlation_form_prior("hierarchical_beta_conjugate"), 4);
  expect_equal(correlation_form_prior("fdhfdhgfdhgf"), NULL);
})

test_that("STAN correlation priors are generated", {
  lkj_param <- 10
  inv_wish_params <- list(2, diag(9))
  beta_params <- list(matrix(5,9,9), matrix(2,9,9))

  STAN_data <-list("TEST_lkj" = array(lkj_param,1),
                   "TEST_wish_nu" = numeric(0),
                   "TEST_wish_sigma" =  matrix(0,0,0),
                   "TEST_beta_1" =  matrix(0,0,0),
                   "TEST_beta_2" =  matrix(0,0,0))
  expect_equal(correlation_prior(params = lkj_param,form=0,type="TEST", allow_hierarchical = FALSE),STAN_data)

  STAN_data[[1]] <- numeric(0)
  STAN_data[[2]] <- array(inv_wish_params[[1]],1)
  STAN_data[[3]] <- inv_wish_params[[2]]

  expect_equal(correlation_prior(params = inv_wish_params,form=1,type="TEST", allow_hierarchical = FALSE),STAN_data)

  STAN_data[[2]] <- numeric(0)
  STAN_data[[3]] <- matrix(0,0,0)
  STAN_data[[4]] <- beta_params[[1]]
  STAN_data[[5]] <- beta_params[[2]]

  expect_equal(correlation_prior(params = beta_params,form=2,type="TEST", allow_hierarchical = FALSE),STAN_data)
})

test_that("STAN correlation priors are generated for hierachical types", {
  #Hierarchical priors
  hierarchical_params <- list(0.8, 0.6, 0.2, 0.7)
  hierarchical_beta_conjugate_params <- list(0.8, 0.8, 0.5)
  STAN_data <-list("TEST_lkj" = numeric(0),
                   "TEST_wish_nu" = numeric(0),
                   "TEST_wish_sigma" =  matrix(0,0,0),
                   "TEST_beta_1" =  matrix(0,0,0),
                   "TEST_beta_2" =  matrix(0,0,0),
                   "TEST_hierarchical_beta_hyper_params" = unlist(hierarchical_params))

  expect_equal(correlation_prior(params = hierarchical_params,form=3,type="TEST", allow_hierarchical = TRUE),STAN_data)

  STAN_data[[6]] <- unlist(hierarchical_beta_conjugate_params)

  expect_equal(correlation_prior(params = hierarchical_beta_conjugate_params,form=4,type="TEST", allow_hierarchical = TRUE),STAN_data)
})

