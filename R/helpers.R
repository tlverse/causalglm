# Helpers

family_list <- list(OR = "binomial", RR = "poisson", CATE = "gaussian")


check_arguments <- function(formula, data, W, A, Y) {
  tryCatch(
    {
      formula <- as.formula(formula)
    },
    error = function(...) {
      stop("Unable to cast `formula` into an R formula object. This should be a character string specifying a valid formula or a formula object.")
    }
  )

  tryCatch(
    {
      data <- as.data.table(data)
    },
    error = function(...) {
      stop("Unable to cast `data` into data.table. Make sure this is a matrix or data.frame.")
    }
  )
  if (!is.character(W)) {
    stop("`W` should be a character vector of baseline covariates.")
  } else if (!(all(W %in% colnames(data)))) {
    stop("Not all variables in `W` were found in `data`.")
  }
  if (length(A) != 1) {
    stop("`A` should be a single character specifying the treatment variable name in `data`.")
  } else if (!(A %in% colnames(data))) {
    stop("Variable `A` was not found in `data`.")
  }
  if (length(Y) != 1) {
    stop("`Y` should be a single character specifying the treatment variable name in `data`.")
  } else if (!(Y %in% colnames(data))) {
    stop("Variable `Y` was not found in `data`.")
  }
}
