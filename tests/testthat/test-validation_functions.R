test_that("Validation for initial input data works", {
  error_message <- paste0("TEST data is not in the correct form. This data should be passed through as a list, the first element being a matrix or data frame containing the data, the second element being the covariance matrix as a matrix or data frame.")
  expect_error(validate_input_data("TEST", NULL), error_message)
  expect_error(validate_input_data("TEST", list(x=1)), error_message)
  expect_error(validate_input_data("TEST", list(x=1, y=2, z=3, t =4)), error_message)
  expect_error(validate_input_data("TEST", list(x=1, y=2)), error_message)
  expect_error(validate_input_data("TEST", list(x=as.matrix(1), y=2)), error_message)
  expect_error(validate_input_data("TEST", list(x=1, y=as.matrix(2))), error_message)
  expect_null(validate_input_data("TEST", list(x=as.matrix(1), y=as.matrix(2))))
  expect_null(validate_input_data("TEST", list(x=as.matrix(1), y=as.matrix(2), "TEST")))
})


test_that("Validation for data/covariance compatibility works", {

  data_matrix <- matrix(0, 40, 9)
  covariance_matrix <- matrix(0,8,9)
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               "TEST covariance matrix is not a positive definite square matrix.")

  covariance_matrix <- matrix(0,8,8)
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               "TEST covariance matrix is not a positive definite square matrix.")

  covariance_matrix <- diag(8)
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               "TEST data has 9 columns but the covariance matrix has 8 columns. Variables for the covariance matrix should match the observed data.")

  covariance_matrix <- diag(9)
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               paste0("TEST data frame variables are not named. Each column of the data frame",
                     " and the covariance matrix should be named appropriately and match each other."))

  colnames(data_matrix) <- c("A","B","C","D","E","F","G","H","I")
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               paste0("TEST covariance matrix variables are not named. Each column of the data frame",
                      " and the covariance matrix should be named appropriately and match each other."))

  colnames(covariance_matrix) <- paste(1:9)
  expect_error(validate_compatibility("TEST", list(data_matrix, covariance_matrix)),
               paste0("TEST covariance matrix variables do not match the associated data frame. Each ",
                      "column of the data frame and covariance matrix should match each other and be ",
                      "in the same order. Data variables: A, B, C, D, E, F, G, H, I Covariance matrix",
                      " variables:  1, 2, 3, 4, 5, 6, 7, 8, 9"))

  colnames(covariance_matrix) <- colnames(data_matrix)
  expect_null(validate_compatibility("TEST", list(data_matrix, covariance_matrix)))
})


test_that("Misc validation for correlation priors works", {

  expect_error(validate_correlation_priors("TEST", "FAKECORRELATIONS", numeric(0), 9),
               paste0("Invalid prior choice on TEST correlation matrix. Correlation priors should be",
                      " one of 'lkj', 'inv_wishart', or 'beta'. Prior choice: fakecorrelations"))
})

test_that("Validation for LKJ correlation priors works", {
  lkj_error <- paste0("Invalid lkj parameter for TEST correlation matrix prior. This should be a ",
                      "list with exactly one real number.")
  expect_error(validate_correlation_priors("TEST", "lKj", numeric(0), 9),lkj_error)
  expect_error(validate_correlation_priors("TEST", "lKj", list(1,2), 9),lkj_error)
  expect_null(validate_correlation_priors("TEST", "lKj", list(1), 9))
})

test_that("Validation for inverse wishart correlation priors works", {
  inv_wishart_error <- "Invalid inverse Wishart parameters for TEST correlation matrix prior."
  expect_error(validate_correlation_priors("TEST", "inv_wishart", numeric(0), 9),inv_wishart_error)
  expect_error(validate_correlation_priors("TEST", "inv_wishart", list(1,2,4), 9),inv_wishart_error)
  expect_error(validate_correlation_priors("TEST", "inv_wishart", list(1), 9),inv_wishart_error)
  expect_error(validate_correlation_priors("TEST", "inv_wishart", list(1,2), 9),inv_wishart_error)

  expect_null(validate_correlation_priors("TEST", "inv_wishart", list(1,diag(9)), 9))
})

test_that("Validation for beta correlation priors works", {
  beta_error <- paste0("Invalid beta parameters for TEST correlation matrix priors. These should be",
                       " square, symmetric matrices with the same dimension as the number of model outputs.")
  expect_error(validate_correlation_priors("TEST", "bEtA", numeric(0), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(1), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(1,2,3), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(1, 2), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(1:2,1,2), 2), 9),beta_error) # Checking matrices
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(1:2,1,2), matrix(1:2,1,2)), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(1:4,2,2), matrix(1:2,1,2)), 9),beta_error)# Checking square matrices
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(1:4,2,2), matrix(1:4,2,2)), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(0,2,2), matrix(1:4,2,2)), 9),beta_error) # Checking symmetric
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(0,2,2), matrix(0,2,2)), 9),beta_error)
  expect_error(validate_correlation_priors("TEST", "bEtA", list(matrix(0,9,9), matrix(0,2,2)), 9),beta_error) # Right dimensionality

  expect_null(validate_correlation_priors("TEST", "bEtA", list(matrix(0,9,9), matrix(0,9, 9)), 9))
})

