


#' contglm
#' Robust generalized linear models for interpretable causal inference for continuous or ordered treatments.
#' Currently, only supports models for the CATE, the log OR, and the log RR of the form: `1(A>0) * f(W) + A * g(W)` with `f` and `g` user-specified.
#' @param formula_continuous An R formula object specifying the continuous component of the parametric form of the continuous treatment CATE.
#' That is (using CATE as example), \code{formula_binary} specifies the interaction with `A` in the model `E[Y|A=a,W] - E[Y|A=0,W] = 1(A>0) * f(W) + A * g(W)`.
#' @param formula_binary An R formula object specifying the binary component of the parametric form of the continuous treatment estimand.
#' That is (using CATE as example), \code{formula_binary} specifies the interaction with `1(A>0)` in the model `E[Y|A=a,W] - E[Y|A=0,W] = 1(A>0) * f(W) + A * g(W)`.
#' By default, the same as \code{formula_continuous}
#' @param data A data.frame or matrix containing the numeric values corresponding with the nodes \code{W}, \code{A} and \code{Y}.
#' Can also be a \code{npglm} fit/output object in which case machine-learning fits are reused (see vignette).
#' @param W A character vector of covariates contained in \code{data}
#' @param A A character name for the treatment assignment variable contained in \code{data}
#' @param Y A character name for the outcome variable contained in \code{data} (outcome can be continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "SuperLearner": A stacked ensemble of all of the below that utilizes cross-validation to adaptivelly choose the best learner.
#' "HAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001}
#' "glm": Fit nuisances with parametric model.
#' "glmnet": Learn using lasso with glmnet.
#' "gam": Learn using generalized additive models with mgcv.
#' "mars": Multivariate adaptive regression splines with \code{earth}.
#' "ranger": Robust random-forests with the package \code{Ranger}
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7.
#' Note speed can vary  depending on learner choice!
#' @param estimand Estimand/parameter to estimate. Options are:
#' `CATE`: conditional treatment effect using working model `CATE(a,W) = E[Y|A=a,W] - E[Y|A=0,W] = 1(a>0) * f(W) + a * g(W)`
#' `OR`: conditional odds ratio using working model `log OR(a,W) = log P(Y=1|A=a,W)/P(Y=0|A=a,W) - log P(Y=1|A=0,W)/P(Y=0|A=0,W)  = 1(a>0) * f(W) + a * g(W)`
#' `RR`: conditional relative risk using working model `log RR(a,W) = log E[Y|A=a,W] - log E[Y|A=0,W] = 1(a>0) * f(W) + a * g(W)`
#' @param cross_fit Whether to cross-fit the initial estimator. This is always set to FALSE if argument \code{sl3_Learner_A} and/or \code{sl3_Learner_Y} is provided.
#' learning_method = `SuperLearner` is always cross-fitted (default).
#'  learning_method = `xgboost` and `ranger` are always cross-fitted regardless of the value of \code{cross_fit}
#'  All other learning_methods are only cross-fitted if `cross_fit=TRUE`.
#'  Note, it is not necessary to cross-fit glm, glmnet, gam or mars as long as the dimension of W is not too high.
#'  In smaller samples and lower dimensions, it may in fact hurt to cross-fit.
#' @param sl3_Learner_A A \code{sl3} Learner object to use to estimate nuisance functions `P(A>0|W)` and `E[A|W]`` with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided.
#' If you wish to cross-fit the learner \code{sl3_Learner_A} then do: sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param sl3_Learner_Y A \code{sl3} Learner object to use to nonparametrically E[Y|A,W] with machine-learning.
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
#'
#'
#' @export
contglm <- function(formula_continuous, formula_binary = formula_continuous, data, W, A, Y, estimand = c("CATE", "OR", "RR"), learning_method = c("HAL", "SuperLearner", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"), cross_fit = FALSE, sl3_Learner_A = NULL, sl3_Learner_Y = NULL, formula_Y = as.formula(paste0("~ . + . *", A)), formula_HAL_Y = paste0("~ . + h(.,", A, ")"), HAL_args_Y = list(smoothness_orders = 1, max_degree = 2, num_knots = c(15, 10, 1)), HAL_fit_control = list(parallel = F), delta_epsilon = 0.025, verbose = TRUE, ...) {
  formula <- NULL

  if (inherits(data, "contglm")) {
    if (verbose) {
      print("Reusing previous fit...")
    }
    args <- data$args
    args$formula <- formula
    tmle3_input <- data$tmle3_input

    tmle3_fit <- refit_glm(data, formula, estimand = estimand, verbose = verbose)
    data <- args$data
  } else {
    check_arguments(formula_continuous, data, W, A, Y)
    check_arguments(formula_binary, data, W, A, Y)
    args <- list(formula_binary = formula_binary, formula_continuous = formula_continuous, formula = formula, data = data, W = W, A = A, Y = Y)


    weights <- NULL

    estimand <- match.arg(estimand)
    learning_method <- match.arg(learning_method)

    if (!is.null(weights)) {
      data$weights <- weights
    } else {
      data$weights <- 1
    }

    superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
      Lrnr_glmnet$new(), Lrnr_glm$new(), Lrnr_gam$new(), Lrnr_earth$new(),
      Lrnr_ranger$new(), Lrnr_xgboost$new(verbose = 0, max_depth = 3), Lrnr_xgboost$new(verbose = 0, max_depth = 4), Lrnr_xgboost$new(verbose = 0, max_depth = 5)
    ), full_fit = T), Lrnr_cv_selector$new(loss_squared_error))

    learner_list_A <- list(
      HAL = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 1, num_knots = c(10, 3)), SuperLearner = superlearner_default, glmnet = Lrnr_glmnet$new(), glm = Lrnr_glm$new(), gam = Lrnr_gam$new(), mars = Lrnr_earth$new(),
      ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new(Stack$new(Lrnr_xgboost$new(verbose = 0, max_depth = 3, eval_metric = "logloss"), Lrnr_xgboost$new(verbose = 0, max_depth = 4, eval_metric = "logloss"), Lrnr_xgboost$new(verbose = 0, max_depth = 5, eval_metric = "logloss")), full_fit = TRUE), Lrnr_cv_selector$new(loss_loglik_binomial))
    )

    if (is.null(sl3_Learner_A)) {
      sl3_Learner_A <- learner_list_A[[learning_method]]
      if (learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
        sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A, full_fit = TRUE)
      }
    }
    binary <- all(Y %in% c(0, 1))
    if (is.null(sl3_Learner_Y)) {
      if (learning_method == "HAL") {
        wrap_in_Lrnr_glm_sp <- FALSE
        # Allow for formula_HAL
        sl3_Learner_Y <- Lrnr_hal9001$new(
          formula_HAL = formula_HAL_Y, family = family_list[[estimand]],
          smoothness_orders = HAL_args_Y$smoothness_orders,
          max_degree = HAL_args_Y$max_degree,
          num_knots = HAL_args_Y$num_knots, fit_control = HAL_fit_control
        )
      } else if (estimand == "RR" && !binary) {
        superlearner_RR <- make_learner(Pipeline, Lrnr_cv$new(list(
          Lrnr_glmnet$new(family = "poisson", formula = formula_Y), Lrnr_glm$new(family = poisson(), formula = formula_Y), Lrnr_gam$new(family = poisson()),
          Lrnr_xgboost$new(verbose = 0, max_depth = 3, objective = "count:poisson"), Lrnr_xgboost$new(verbose = 0, max_depth = 4, objective = "count:poisson"), Lrnr_xgboost$new(verbose = 0, max_depth = 5, objective = "count:poisson")
        ), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))

        learner_list_Y0W_RR <- list(
          SuperLearner = superlearner_RR, glmnet = Lrnr_glmnet$new(formula = formula_Y, family = "poisson"), glm = Lrnr_glm$new(formula = formula_Y, family = poisson()), gam = Lrnr_gam$new(family = poisson()),
          xgboost = make_learner(Pipeline, Lrnr_cv$new(Stack$new(Lrnr_xgboost$new(verbose = 0, max_depth = 3, objective = "count:poisson", eval_metric = "error"), Lrnr_xgboost$new(verbose = 0, max_depth = 4, objective = "count:poisson", eval_metric = "error"), Lrnr_xgboost$new(verbose = 0, max_depth = 5, objective = "count:poisson", eval_metric = "error")), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
        )
        sl3_Learner_Y <- learner_list_Y0W_RR[[learning_method]]
      } else {
        superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(list(
          Lrnr_glmnet$new(formula = formula_Y), Lrnr_glm$new(formula = formula_Y), Lrnr_gam$new(), Lrnr_earth$new(formula = formula_Y),
          Lrnr_ranger$new(), Lrnr_xgboost$new(verbose = 0, max_depth = 3), Lrnr_xgboost$new(verbose = 0, max_depth = 4), Lrnr_xgboost$new(verbose = 0, verbose = 0, max_depth = 5)
        ), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
        learner_list_Y <- list(
          SuperLearner = superlearner_default, glmnet = Lrnr_glmnet$new(formula = formula_Y), glm = Lrnr_glm$new(formula = formula_Y), gam = Lrnr_gam$new(), mars = Lrnr_earth$new(formula = formula_Y),
          ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new(Stack$new(Lrnr_xgboost$new(verbose = 0, verbose = 0, max_depth = 3), Lrnr_xgboost$new(verbose = 0, verbose = 0, max_depth = 4), Lrnr_xgboost$new(verbose = 0, max_depth = 5)), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
        )
        sl3_Learner_Y <- learner_list_Y[[learning_method]]
      }
      if (learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
        sl3_Learner_Y <- Lrnr_cv$new(sl3_Learner_Y, full_fit = TRUE)
      }
    }

    if (estimand == "CATE") {
      fA <- function(x) {
        x
      }
      dfA <- function(x) {
        1
      }
      submodel <- NULL
    } else if (estimand == "RR") {
      fA <- function(x) {
        log(x)
      }
      dfA <- function(x) {
        1 / x
      }
      submodel <- NULL
    } else if (estimand == "OR") {
      fA <- function(x) {
        log(x) - log(1 - x)
      }
      dfA <- function(x) {
        1 / (x * (1 - x))
      }
      submodel <- "binomial"
    }

    tmle_spec_np <- tmle3_Spec_contglm$new(
      formula_continuous = formula_continuous, formula_binary = formula_binary,
      fA = fA, dfA = dfA,
      submodel = submodel,
      delta_epsilon = delta_epsilon, verbose = verbose, include_A_binary = TRUE
    )

    # tmle_spec_np <- tmle3_Spec_contCATE$new(formula_continuous = formula_continuous, formula_binary = formula_binary, delta_epsilon = delta_epsilon, verbose = verbose, include_A_binary = TRUE)
    learner_list <- list(A = sl3_Learner_A, A_binary = sl3_Learner_A, Y = sl3_Learner_Y)
    node_list <- list(W = W, A = A, Y = Y)

    tmle3_input <- list(tmle_spec = tmle_spec_np, data = data, node_list = node_list, learner_list = learner_list, delta_epsilon = delta_epsilon)
    tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec_np, data, node_list, learner_list)))
  }
  coefs <- tmle3_fit$summary
  coefs <- coefs[, -3]
  if (estimand %in% c("CATE", "CATT", "TSM")) {
    coefs <- coefs[, 1:6]
  } else {
    cur_names <- colnames(coefs)
    cur_names <- gsub("transformed", "exp", cur_names)
    colnames(coefs) <- cur_names
  }
  n <- nrow(data)
  Zscore <- abs(sqrt(n) * coefs$tmle_est / coefs$se)
  pvalue <- signif(2 * (1 - pnorm(Zscore)), 5)
  coefs$Z_score <- Zscore
  coefs$p_value <- pvalue

  tmp <- coefs$param
  if (estimand %in% c("OR", "RR")) {
    formula_fit <- paste0("log ", coefs$type[1], "(W) = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
  } else {
    formula_fit <- paste0(coefs$type[1], "(A,W) = ", paste0(signif(coefs$tmle_est, 3), " * ", tmp, collapse = " + "))
  }

  output <- list(estimand = estimand, formula_fit = formula_fit, coefs = coefs, tmle3_fit = tmle3_fit, tmle3_input = tmle3_input, args = args)
  class(output) <- c("contglm", "causalglm")
  return(output)
}
