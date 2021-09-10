
#' causalglmnet
#' High dimensional semiparametric generalized linear models for causal inference using the LASSO.
#' Supports flexible semiparametric conditional average treatment effect (CATE), conditional odds ratio (OR), and conditional relative risk (RR) estimation
#' \code{\link[glmnet]{cv.glmnet}} is used to fit all nuisance parameters. The parametric component of the semiparametric model is not penalized.
#' This function is almost just a wrapper for \code{causalglm}.
#'
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param data A data.frame or matrix containing the numeric values corresponding with the nodes \code{W}, \code{A} and \code{Y}.
#' @param W A character vector of covariates contained in \code{data}
#' @param A A character name for the treatment assignment variable contained in \code{data}
#' @param Y A character name for the outcome variable contained in \code{data} (outcome can be continuous, nonnegative or binary depending on method)
#' @param estimand Estimand/parameter to estimate. Choices are:
#' CATE: Estimate conditional average treatment effect with \code{\link[tmle3]{Param_spCATE}} assuming it satisfies parametric model \code{formula}.
#' OR: Estimate conditional odds ratio with \code{\link[tmle3]{Param_spOR}} assuming it satisfies parametric model \code{formula}.
#' RR: Estimate conditional relative risk with \code{\link[tmle3]{Param_spRR}} assuming it satisfies parametric model \code{formula}.
#' @param cross_fit Whether to cross-fit the initial estimator. This is always set to FALSE if argument \code{sl3_Learner} is provided.
#' learning_method = `SuperLearner` is always cross-fitted (default).
#'  learning_method = `xgboost` and `ranger` are always cross-fitted regardless of the value of \code{cross_fit}
#'  All other learning_methods are only cross-fitted if `cross_fit=TRUE`.
#'  Note, it is not necessary to cross-fit glm, glmnet, gam or mars as long as the dimension of W is not too high.
#'  In smaller samples and lower dimensions, it may fact hurt to cross-fit.
#' @param weights An optional vector of weights to use in procedure.
#' @param parallel See \code{\link[glmnet]{cv.glmnet}}
#' @param ... Other arguments to pass to \code{\link[glmnet]{cv.glmnet}}
#' @export
causalglmnet <- function(formula, data, W, A, Y, estimand = c("CATE", "OR", "RR"), max_degree = 1, cross_fit = TRUE, constant_variance_CATE = FALSE, weights = NULL, parallel = TRUE, verbose = TRUE, ...) {
  append_interaction_matrix <- TRUE
  estimand <- match.arg(estimand)

  # HAL_args_Y0W = list(smoothness_orders = 1, max_degree = 1, num_knots = 1)
  HAL_fit_control <- list(parallel = parallel, ...)

  data <- as.data.table(data)
  if (!is.null(weights)) {
    data$weights <- weights
  } else {
    data$weights <- 1
  }
  sl3_Learner_A <- Lrnr_glmnet$new()
  if (constant_variance_CATE) {
    sl3_Learner_var_Y <- Lrnr_mean$new()
  } else {
    sl3_Learner_var_Y <- Lrnr_glmnet$new(formula = formula(paste0("~ . + .*", A)), family = "poisson")
  }
  sl3_Learner_Y <- Lrnr_hal9001_semiparametric$new(
    formula_sp = formula, family = family_list[[estimand]],
    interaction_variable = A,
    smoothness_orders = 1,
    max_degree = max_degree,
    num_knots = 1, fit_control = HAL_fit_control
  )

  if (cross_fit) {
    sl3_Learner_Y <- Lrnr_cv$new(sl3_Learner_Y, full_fit = TRUE)
    sl3_Learner_var_Y <- Lrnr_cv$new(sl3_Learner_var_Y, full_fit = TRUE)
    sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A, full_fit = TRUE)
  }
  tmle_spec_sp <- tmle3_Spec_spCausalGLM$new(formula = formula, estimand = estimand, append_interaction_matrix = TRUE, wrap_in_Lrnr_glm_sp = FALSE, binary_outcome = FALSE, delta_epsilon = 1, verbose = verbose)
  learner_list <- list(A = sl3_Learner_A, Y = sl3_Learner_Y)
  if (estimand == "CATE") {
    learner_list$var_Y <- sl3_Learner_var_Y
  }
  node_list <- list(weights = "weights", W = W, A = A, Y = Y)
  tmle3_input <- list(tmle_spec_sp = tmle_spec_sp, data = data, node_list = node_list, learner_list = learner_list)
  tmle3_fit <- ((tmle3(tmle_spec_sp, data, node_list, learner_list)))

  output <- list(coefs = tmle3_fit$summary, tmle3_fit = tmle3_fit, tmle3_input = tmle3_input)
  class(output) <- c("causalglmnet", "causalglm")
  return(output)
}
