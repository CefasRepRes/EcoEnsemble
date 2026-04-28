# tests/testthat/test-warning-filter.R

capture_warnings <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    expr,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  msgs
}

test_that("policy: only max-treedepth warnings are suppressed", {
  fake_sampling <- function(...) {
    warning(
      "There were 5 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10.",
      call. = FALSE
    )
    structure(list(sampler = "fake"))
  }

  local_mocked_bindings(
    sampling = fake_sampling,
    .package = "rstan"
  )

  msgs <- capture_warnings({
    fit <- EcoEnsemble:::stan_sampling_with_filter(
      mod = NULL,
      data = list(),
      control = list()
    )
    expect_type(fit, "list")
  })

  expect_length(msgs, 0)
})

test_that("policy: mixed warnings (treedepth + others) re-emit all, in order, deduped", {
  fake_sampling <- function(...) {
    warning(
      "There were 5 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10.",
      call. = FALSE
    )
    warning(
      "Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.",
      call. = FALSE
    )
    warning(
      "Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.",
      call. = FALSE
    )
    warning( # duplicate, should be deduped
      "There were 5 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10.",
      call. = FALSE
    )
    structure(list(sampler = "fake"))
  }

  local_mocked_bindings(
    sampling = fake_sampling,
    .package = "rstan"
  )

  msgs <- capture_warnings({
    fit <- EcoEnsemble:::stan_sampling_with_filter(
      mod = NULL,
      data = list(),
      control = list()
    )
    expect_type(fit, "list")
  })

  expect_gte(length(msgs), 3L)
  expect_equal(sum(grepl("maximum treedepth", msgs, ignore.case = TRUE)), 1L)
  expect_true(any(grepl("maximum treedepth", msgs[1], ignore.case = TRUE)))
  expect_true(any(grepl("^Bulk Effective Samples Size", msgs[2])))
  expect_true(any(grepl("^Tail Effective Samples Size", msgs[3])))
})

test_that("policy: non-treedepth warnings are emitted", {
  fake_sampling <- function(...) {
    warning("The largest R-hat is NA, indicating chains have not mixed.", call. = FALSE)
    structure(list(sampler = "fake"))
  }

  local_mocked_bindings(
    sampling = fake_sampling,
    .package = "rstan"
  )

  msgs <- capture_warnings({
    fit <- EcoEnsemble:::stan_sampling_with_filter(
      mod = NULL,
      data = list(),
      control = list()
    )
    expect_type(fit, "list")
  })

  expect_length(msgs, 1)
  expect_match(msgs[[1]], "^The largest R-hat is NA")
})

