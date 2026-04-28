stan_sampling_with_filter <- function(mod, data, control = NULL, ...) {
  warn_buf <- character(0)

  fit <- withCallingHandlers(
    rstan::sampling(mod, data = data, control = control, ...),
    warning = function(w) {
      warn_buf <<- c(warn_buf, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (length(warn_buf) > 0) {
    # preserve order, drop duplicates
    msgs <- warn_buf[!duplicated(warn_buf)]
    emit_msgs <- msgs
    # robust match for Stan treedepth warnings
    is_treedepth <- grepl(
      "exceeded the maximum treedepth|maximum\\s+treedepth(.*exceeded)?",
      msgs, ignore.case = TRUE, perl = TRUE
    )
    is_pairs <- grepl(
      "examine\\s+the\\s+pairs\\(\\)\\s+plot",
      msgs, ignore.case = TRUE, perl = TRUE
    )
    is_ess <- grepl(
      "effective\\s+sample\\s+size|\\bESS\\b.*too\\s+low|bulk\\s+effective\\s+sample\\s+size|tail\\s+effective\\s+sample\\s+size",
      msgs, ignore.case = TRUE, perl = TRUE
    )
    emit_msgs[is_ess] <- paste0(
      emit_msgs[is_ess],
      ' Run - vignette("ESS_vignette", package = "EcoEnsemble").'
    )

    suppressible <- is_treedepth | is_pairs

    if (!all(suppressible)) {
      for (msg in emit_msgs) {
        warning(msg, call. = FALSE)
      }
    }
  }

  return(fit)
}
