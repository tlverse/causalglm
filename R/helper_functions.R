#' Plotting of univariate marginal structural models
#' @param object A `msmglm` object
#' @param data The data containing `V` for plotting
#' @import ggplot2
#' @export
plot_msm <- function(object, data = object$args$data) {
  require(ggplot2)
  if (!inherits(object, "msmglm")) {
    stop("object must be a msmglm object.")
  }
  if (length(object$args$V) > 1) {
    stop("MSM plotting is only supported for one-dimensional `V`.")
  }
  preds_full <- predict(object, data = data)
  index <- which(colnames(preds_full) == "se") - 1
  x <- data[[object$args$V]]
  y <- preds_full[[index]]
  data <- data.frame(x = x, y = y, upper = preds_full$CI_right, lower = preds_full$CI_left)

  plot <- ggplot(data, aes(x = x, y = y)) +
    geom_smooth(se = F) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    xlab(paste0("V = ", object$args$V)) +
    ylab("MSM(V)") +
    ggtitle(object$formula_fit) +
    theme(plot.title = element_text(size = 8))
  return(plot)
}



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
  out <- object$coefs
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
  preds_new <- as.data.frame(preds_new)
  return(preds_new)
}


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
