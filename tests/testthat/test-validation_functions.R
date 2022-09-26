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

