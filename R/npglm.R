


#' npglm
#' Nonparametric robust generalized linear models for causal inference and marginal structural models (for some estimands).
#' Supports flexible nonparametric conditional average treatment effect (CATE), conditional odds ratio (OR), and conditional relative risk (RR) estimation,
#' ... where a user-specified working parametric model for the estimand is viewed as an approximation of the true estimand and nonparametrically correct inference is given for these approximations.
#' Specifically, a causal projection of the true estimand onto the working-model is estimated; the parametric model is not assumed correct.
#' The estimates and inference obtained by `npglm` are robust and nonparametrically correct, which comes at a small cost in confidence interval width relative to `spglm`.
#' Highly Adaptive Lasso (HAL) (see \code{\link[hal9001]{fit_hal}}), a flexible and adaptive spline regression estimator, is recommended for medium-small to large sample sizes.
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param data A data.frame or matrix containing the numeric values corresponding with the nodes \code{W}, \code{A} and \code{Y}.
#' @param W A character vector of covariates contained in \code{data}
#' @param A A character name for the treatment assignment variable contained in \code{data}
#' @param Y A character name for the outcome variable contained in \code{data} (outcome can be continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "SuperLearner": A stacked ensemble of all of the below that utilizes cross-validation to adaptivelly choose the best learner.
#' "HAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001} Good for most sample sizes when propertly tuned. See arguments \code{max_degree_Y0W} and \code{num_knots_Y0W}.
#' "glm": Fit nuisances with parametric model. Best for smaller sample sizes (e.g. n =30-100). See arguments \code{glm_formula_A}, \code{glm_formula_Y} and \code{glm_formula_Y0}.
#' "glmnet": Learn using lasso with glmnet. Best for smaller sample sizes (e.g. n =30-100)
#' "gam": Learn using generalized additive models with mgcv. Good for small-to-medium-small sample sizes.
#' "mars": Multivariate adaptive regression splines with \code{earth}. Good for small-to-medium-small sample sizes.
#' "ranger": Robust random-forests with the package \code{Ranger} Good for medium-to-large sample sizes.
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7. Good for medium-to-large sample sizes.
#' We recommend performing simulations checking 95% CI coverage when choosing learners (especially in smaller sample sizes).
#' Note speed can vary significantly depending on learner choice!
#' @param estimand Estimand/parameter to estimate. Choices are:
#' `CATE`: Estimate the best parametric approximation of the conditional average treatment effect with \code{\link[tmle3]{Param_npCATE}} using the parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{npglm} returns a nonparametric and efficient estimator for the ATE (Marginal average treatment effect).
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model.
#' `CATT`: Estimate the best parametric approximation of the conditional average treatment effect among the treated with \code{\link[tmle3]{Param_npCATE}} using the parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{npglm} returns a nonparametric and efficient estimator for the ATT (Marginal average treatment effect among the treated).
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model using only the observations with `A=1` (among the treated).
#' `TSM`: Estimate the best parametric approximation of the conditional treatment-specific mean `E[Y|A=a,W]` for `a` in \code{levels_A}.
#' Note: if \code{formula} = `~1` then \code{npglm} returns a nonparametric and efficient estimator for the TSM (Marginal treatment-specific mean).
#' Specifically, this estimand is the least-squares projection of the true TSM onto the parametric working model.
#' `OR`: Estimate the best parametric approximation of the conditional odds ratio with \code{\link[tmle3]{Param_npOR}} using the parametric model \code{formula}.
#' Specifically, this estimand is the log-likelihood projection of the true conditional odds ratio onto the partially-linear logistic regression model with the true `E[Y|A=0,W]` used as offset.
#' `RR`: Projection of the true conditional relative risk onto a exponential working-model using log-linear/poisson regression.
#' Note: if \code{formula} = `~1` then \code{npglm} returns a nonparametric and efficient estimator for the marginal relative risk (`E_W[E[Y|A=1,W]]/E_W[E[Y|A=0,W]]``).
#' @param levels_A Only used if \code{estimand} = `TSM` in which case the TSM is learned for all levels in \code{levels_A}.
#' @param cross_fit Whether to cross-fit the initial estimator. This is always set to FALSE if argument \code{sl3_Learner} is provided.
#' learning_method = `SuperLearner` is always cross-fitted (default).
#'  learning_method = `xgboost` and `ranger` are always cross-fitted regardless of the value of \code{cross_fit}
#'  All other learning_methods are only cross-fitted if `cross_fit=TRUE`.
#'  Note, it is not necessary to cross-fit glm, glmnet, gam or mars as long as the dimension of W is not too high.
#'  In smaller samples and lower dimensions, it may fact hurt to cross-fit.
#' @param sl3_Learner_A A \code{sl3} Learner object to use to estimate nuisance function P(A=1|W) with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided.
#' If you wish to cross-fit the learner \code{sl3_Learner} then do: sl3_Learner <- Lrnr_cv$new(sl3_Learner).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param sl3_Learner_Y A \code{sl3} Learner object to use to nonparametrically [Y|A,W] with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided.
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param formula_Y Only used if `learning_method %in% c("glm", "earth", "gam", "glmnet")`. A R \code{formula} object that specifies the design matrix to be passed to the Learner specified by learning_method: "glm", "earth", "gam", "glmnet".
#' By default, `formula_Y = . + A*.` so that additive learners still model treatment interactions.
#' @param formula_HAL_Y A HAL formula string to be passed to \code{\link[hal9001]{fit_hal}}). See the `formula` argument of \code{\link[hal9001]{fit_hal}}) for syntax and example use.
#' @param weights An optional vector of weights to use in procedure.
#' @param HAL_args_Y A list of parameters for the semiparametric Highly Adaptive Lasso estimator for E[Y|A,W].
#' Possible parameters are:
#' 1. `smoothness_orders`: Smoothness order for HAL estimator of E[Y|A,W] (see \code{\link[hal9001]{fit_hal}})
#' smoothness_order_Y0W = 1 is piece-wise linear. smoothness_order_Y0W = 0 is piece-wise constant.
#' 2. `max_degree`: Max interaction degree for HAL estimator of E[Y|A,W] (see \code{\link[hal9001]{fit_hal}})
#' 3. `num_knots`: A vector of the number of knots by interaction degree for HAL estimator of E[Y|A=0,W] (see \code{\link[hal9001]{fit_hal}}). Used to generate spline basis functions.
#' num_knots_Y0W = c(1) is equivalent to main term glmnet/LASSO. (Assuming max_degree_Y0W = 1)
#' num_knots_Y0W = c(1,1) is equivalent to glmnet/LASSO with both main-terms and all two-way interactions (e.g. Y~ W1 + A + A W1+ W1 + W1*W2 + ...).  (Assuming max_degree_Y0W = 2)
#' num_knots_Y0W = c(10) is an additive piece-wise linear model with 10 knot points.  (Assuming max_degree_Y0W = 1)
#' num_knots_Y0W = c(10,5) is a bi-additive model with the same one-way basis functions as above, but also two-way interaction piece-wise linear basis functions generated by main-term one-way basis functions with 5 knots. (Assuming max_degree_Y0W = 2)
#' num_knots_Y0W = c(10,5,1) generates same basis functions as above and also the three-way interaction basis functions with only a single knot point at the origin (e.g. the triple interaction `W1*W2*W3`) (Assuming max_degree_Y0W = 3)
#' @param HAL_fit_control See the argument `fit_control` of (see \code{\link[hal9001]{fit_hal}}).
#' @param delta_epsilon Step size of iterative targeted maximum likelihood estimator. `delta_epsilon = 1 ` leads to large step sizes and fast convergence. `delta_epsilon = 0.005` leads to slower convergence but possibly better performance.
#' Useful to set to a large value in high dimensions.
#' @param ... Other arguments to pass to main routine (spCATE, spOR, spRR)
#' @export
npglm <- function(formula, data, W, A, Y, estimand = c("CATE", "CATT", "TSM", "OR", "RR"), learning_method = c("HAL", "SuperLearner", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"), levels_A = sort(unique(data[[A]])), cross_fit = FALSE, sl3_Learner_A = NULL, sl3_Learner_Y = NULL , formula_Y = as.formula(paste0("~ . + . *", A)), formula_HAL_Y = paste0("~ . + h(.,", A, ")"), HAL_args_Y = list(smoothness_orders = 1, max_degree = 2, num_knots = c(15, 10, 1)), HAL_fit_control = list(parallel = F), delta_epsilon = 0.025, verbose = TRUE, ...) {
  check_arguments(formula, data, W, A, Y)
  tryCatch({
  data <- as.data.table(data)
  }, error = function(...){
    stop("Unable to cast `data` into data.table. Make sure this is a matrix or data.frame.")
  })
  if(!is.character(W)) {
    stop("`W` should be a character vector of baseline covariates.")
  } else if (!(all(W %in% colnames(data)))) {
    stop("Not all variables in `W` were found in `data`.")
  }
  if(length(A)!=1){
    stop("`A` should be a single character specifying the treatment variable name in `data`.")
  } else if (!(A %in% colnames(data))) {
    stop("Variable `A` was not found in `data`.")
  }
  if(length(Y)!=1){
    stop("`Y` should be a single character specifying the treatment variable name in `data`.")
  } else if (!(Y %in% colnames(data))) {
    stop("Variable `Y` was not found in `data`.")
  }
  weights = NULL
  
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
  if (estimand != "TSM") {
    levels_A <- max(as.numeric(levels_A))
  }
  tmle_spec_np <- tmle3_Spec_npCausalGLM$new(formula = formula, estimand = estimand, delta_epsilon = delta_epsilon, verbose = verbose, treatment_level = levels_A)
  learner_list <- list(A = sl3_Learner_A, Y = sl3_Learner_Y)
  node_list <- list(weights = "weights", W = W, A = A, Y = Y)

  tmle3_input <- list(tmle_spec_np = tmle_spec_np, data = data, node_list = node_list, learner_list = learner_list)
  tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec_np, data, node_list, learner_list)))

  output <- list(coefs = tmle3_fit$summary, tmle3_fit = tmle3_fit, tmle3_input = tmle3_input)
  class(output) <- c("npglm", "causalglm")
  return(output)
}
