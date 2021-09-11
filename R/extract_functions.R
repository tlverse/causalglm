#' @export
summary.causalglm <- function(object) {
  print(object)
  cat("\n\n")
  cat("Coefficient estimates and inference:")
  cat("\n")
  print(object$coefs)
  return(invisible(object$coefs))
}

#' @export
coef.causalglm <- function(object) {
  out <- (object$coefs)
  out
}

#' @export
print.causalglm <- function(object) {
  cat(paste0("A causalglm fit object obtained from ", class(object)[1], " for the estimand ", object$estimand, " with formula: \n"))
  cat(object$formula_fit)
}


#' @export
predict.causalglm <- function(object, data = object$args$data, transformed = TRUE) {
  W <- object$args$W
  formula <- object$args$formula
  estimand <- object$estimand
  V <- model.matrix(formula, as.data.frame(data))
  n <- nrow(object$args$data)

  estimates <- object$coefs$tmle_est
  var_mat <- var(object$tmle3_fit$estimates[[1]]$IC)

  est_grid <- V %*% estimates



  se_grid <- apply(V, 1, function(m) {
    sqrt(sum(m * (var_mat %*% m)))
  })
  Zvalue <- abs(sqrt(n) * est_grid / se_grid)
  pvalue <- signif(2 * (1 - pnorm(Zvalue)), 5)

  if (estimand %in% c("OR", "RR")) {
    linkinv <- exp
  } else {
    linkinv <- function(x) x
  }
  ci <- cbind(est_grid - 1.96 * se_grid / sqrt(n), est_grid + 1.96 * se_grid / sqrt(n))
  if (transformed) {
    ci <- linkinv(ci)
    est_grid <- linkinv(est_grid)
  }

  preds_new <- cbind(V, est_grid, se_grid, ci, Zvalue, pvalue)
  name <- paste0(object$estimand, "(W)")
  if (!transformed && estimand %in% c("OR", "RR")) {
    name <- paste0("log ", name)
  }
  colnames(preds_new) <- c(colnames(V), name, "se", "CI_left", "CI_right", "Z-score", "p-value")
  return(preds_new)
}
