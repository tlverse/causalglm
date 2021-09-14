
#' Marginal structural generalized linear models with robust inference based on nonparametric projections.
#' Nonparametric robust inference for marginal structural models for the CATE, CATT, TSM, and RR.
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param data A data.frame or matrix containing the numeric values corresponding with the nodes \code{W}, \code{A} and \code{Y}.
#' @param V The marginal covariate (or covariate vector) of interest for the marginal structural model.
#' @param W A character vector of covariates/confounders to adjust for contained in \code{data}. \code{V} is automatically added to \code{W}.
#' @param A A character name for the treatment assignment variable contained in \code{data}
#' @param Y A character name for the outcome variable contained in \code{data} (outcome can be continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "SuperLearner": A stacked ensemble of all of the below that utilizes cross-validation to adaptively choose the best learner.
#' "HAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001}. See argument \code{HAL_args_Y}.
#' "glm": Fit nuisances with parametric model.
#' "glmnet": Learn using lasso with glmnet.
#' "gam": Learn using generalized additive models with mgcv.
#' "mars": Multivariate adaptive regression splines with \code{earth}.
#' "ranger": Robust random-forests with the package \code{Ranger}
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7.
#' We recommend performing simulations checking 95% CI coverage when choosing learners (especially in smaller sample sizes).
#' Note speed can vary significantly depending on learner choice!
#' @param estimand Estimand/parameter to estimate. Choices are:
#' `CATE`: Provides nonparametrically-robust inference for the user-specified marginal structural model of `E[CATE(W)|V] := E[E[Y|A=1,W] - E[Y|A=0,W]|V]`.
#' `CATT`: Provides nonparametrically-robust inference for the user-specified marginal structural model of `E[CATE(W)|V, A=1] := E[E[Y|A=1,W] - E[Y|A=0,W]|V, A=1]`.
#' `TSM`: Provides nonparametrically-robust inference for the user-specified marginal structural model of  `E[E[Y|A=a,W]|V]`.
#' `RR`: Provides nonparametrically-robust inference for the user-specified marginal structural relative risk model of  `E[E[Y|A=1,W]|V] / E[E[Y|A=0,W]|V]`.
#' @param levels_A Only used if \code{estimand} = `TSM` in which case the TSM is learned for all levels in \code{levels_A}.
#' @param cross_fit Whether to cross-fit the initial estimator. This is always set to FALSE if argument \code{sl3_Learner_A} and/or \code{sl3_Learner_Y} is provided.
#' learning_method = `SuperLearner` is always cross-fitted (default).
#'  learning_method = `xgboost` and `ranger` are always cross-fitted regardless of the value of \code{cross_fit}
#'  All other learning_methods are only cross-fitted if `cross_fit=TRUE`.
#'  Note, it is not necessary to cross-fit glm, glmnet, gam or mars as long as the dimension of W is not too high.
#'  In smaller samples and lower dimensions, it may in fact hurt to cross-fit.
#' @param sl3_Learner_A A \code{sl3} Learner object to use to estimate nuisance function P(A=1|W) with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided.
#' If you wish to cross-fit the learner \code{sl3_Learner_A} then do: sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param sl3_Learner_Y A \code{sl3} Learner object to use to nonparametrically [Y|A,W] with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided.
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param formula_Y Only used if `learning_method %in% c("glm", "earth", "glmnet")`. A R \code{formula} object that specifies the design matrix to be passed to the Learner specified by learning_method: "glm", "earth", "glmnet".
#' By default, `formula_Y = . + A*.` so that additive learners still model treatment interactions.
#' @param formula_HAL_Y A HAL formula string to be passed to \code{\link[hal9001]{fit_hal}}). See the `formula` argument of \code{\link[hal9001]{fit_hal}}) for syntax and example use.
#' @param HAL_args_Y A list of parameters for the semiparametric Highly Adaptive Lasso estimator for E[Y|A,W].
#' Should contain the parameters:
#' 1. `smoothness_orders`: Smoothness order for HAL estimator of E[Y|A,W] (see \code{\link[hal9001]{fit_hal}})
#' smoothness_order_Y0W = 1 is piece-wise linear. smoothness_order_Y0W = 0 is piece-wise constant.
#' 2. `max_degree`: Max interaction degree for HAL estimator of E[Y|A,W] (see \code{\link[hal9001]{fit_hal}})
#' 3. `num_knots`: A vector of the number of knots by interaction degree for HAL estimator of E[Y|A=0,W] (see \code{\link[hal9001]{fit_hal}}). Used to generate spline basis functions.
#' @param HAL_fit_control See the argument `fit_control` of (see \code{\link[hal9001]{fit_hal}}).
#' @param delta_epsilon Step size of iterative targeted maximum likelihood estimator. `delta_epsilon = 1 ` leads to large step sizes and fast convergence. `delta_epsilon = 0.01` leads to slower convergence but possibly better performance.
#' Useful to set to a large value in high dimensions.
#' @param verbose Passed to \code{tmle3} routines. Prints additional information if TRUE.
#' @param ... Not used
#' @export
msmglm <- function(formula, data, V, W = V, A, Y, estimand = c("CATE", "CATT", "TSM", "RR"), learning_method = c("HAL", "SuperLearner", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"), treatment_level = max(data[[A]]), control_level = min(data[[A]]), cross_fit = FALSE, sl3_Learner_A = NULL, sl3_Learner_Y = NULL, formula_Y = as.formula(paste0("~ .^2")), formula_HAL_Y = NULL, HAL_args_Y = list(smoothness_orders = 1, max_degree = 2, num_knots = c(15, 10, 1)), HAL_fit_control = list(parallel = F), delta_epsilon = 0.025, verbose = TRUE, ...) {
  if (!inherits(data, "msmglm")) {
    if (is.null(formula_HAL_Y)) {
      formula_HAL_Y <- paste0("~ . + h(.,", A, ")")
    }
    if (length(unique(data[[A]])) > 2) {
      formula_HAL_Y <- paste0("~ . + h(.,.)")
    }
    if (missing(V)) {
      V <- data$args$V
    }
    W <- data$args$W
    W <- union(W, V)
    tryCatch(
      {
        formula <- as.formula(formula)
        check <- model.matrix(formula, as.data.frame(data)[, V, drop = F])
      },
      error = function(...) {
        print(...)
        stop("`formula` should specify a marginal structural model and only depend on variables in `V`.")
      }
    )
  }

  output <- npglm(formula, data, W, A, Y, estimand, learning_method, treatment_level, control_level, cross_fit, sl3_Learner_A, sl3_Learner_Y, formula_Y, formula_HAL_Y, HAL_args_Y, HAL_fit_control, delta_epsilon, verbose, ...)

  if (estimand %in% c("TSM")) {
    levels <- output$levels_A
    for (i in seq_along(levels)) {
      item <- output[[i]]
      coefs <- item$coefs
      tmp <- coefs$param
      formula_fit <- paste0("E[TSM(W)|V]", " = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
      item$formula_fit <- formula_fit
      item$args$V <- V
      class(item) <- c("msmglm", "causalglm")
      output[[i]] <- item
    }
    return(output)
  }
  class(output) <- c("msmglm", "causalglm")
  output$args$V <- V
  coefs <- output$coefs
  tmp <- coefs$param
  if (estimand %in% c("OR", "RR")) {
    formula_fit <- paste0("log E[", coefs$type[1], "(W)|V] = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
  } else if (estimand %in% c("CATE")) {
    formula_fit <- paste0("E[CATE(W)|V]", " = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
  } else if (estimand %in% c("CATT")) {
    formula_fit <- paste0("E[CATE(W)|V, A=1]", " = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
  }
  output$formula_fit <- formula_fit

  return(output)
}
