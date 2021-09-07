
#' causalGLMnet
#' High dimensional semiparametric generalized linear models for causal inference using the LASSO.
#' Supports flexible semiparametric conditional average treatment effect (CATE), conditional odds ratio (OR), and conditional relative risk (RR) estimation
#' \code{\link[glmnet]{cv.glmnet}} is used to fit all nuisance parameters. The parametric component of the semiparametric model is not penalized.
#' This function is almost just a wrapper for \code{causalGLM}.
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
causalGLMnet <- function(formula, data, W, A, Y, estimand = c("CATE", "OR", "RR"),    cross_fit = TRUE,  constant_variance_CATE = FALSE, weights = NULL,  parallel = TRUE , verbose = TRUE, ... ){
  append_interaction_matrix = TRUE
  estimand <- match.arg(estimand)
  
  #HAL_args_Y0W = list(smoothness_orders = 1, max_degree = 1, num_knots = 1)
  HAL_fit_control = list(parallel = parallel, ...)
  
  data <- as.data.table(data)
  if(!is.null(weights)) {
    data$weights <- weights
  } else {
    data$weights <- 1
  }
  sl3_Learner_A <- Lrnr_glmnet$new()
  if(constant_variance_CATE) {
    sl3_Learner_var_Y <- Lrnr_mean$new()
  } else {
    sl3_Learner_var_Y <- Lrnr_glmnet$new(formula = formula(paste0("~ . + .*", A)), family = "poisson")
  }
  sl3_Learner_Y <- Lrnr_hal9001_semiparametric$new(formula =  formula, family = family_list[[estimand]], 
                                                       interaction_variable = A,
                                                       smoothness_orders = 1,
                                                       max_degree = 1,
                                                       num_knots = 1, fit_control = HAL_fit_control
      )
   
  if(cross_fit) {
    sl3_Learner_Y <- Lrnr_cv$new(sl3_Learner_Y, full_fit = TRUE)
    sl3_Learner_var_Y <- Lrnr_cv$new(sl3_Learner_var_Y, full_fit = TRUE)
    sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A, full_fit = TRUE)
  }
  tmle_spec_sp <- tmle3_Spec_spCausalGLM$new(formula = formula, estimand = estimand, append_interaction_matrix = TRUE, wrap_in_Lrnr_glm_sp = FALSE , binary_outcome = FALSE, delta_epsilon = 1, verbose = verbose )
  learner_list <- list(A = sl3_Learner_A, Y = sl3_Learner_Y)
  if(estimand == "CATE") {
    learner_list$var_Y <- sl3_Learner_var_Y
  }
  node_list <- list(weights = "weights", W = W, A = A, Y = Y)
  tmle3_input <- list(tmle_spec_sp = tmle_spec_sp, data = data, node_list = node_list,  learner_list = learner_list)
  tmle3_fit <-  ( (tmle3(tmle_spec_sp, data, node_list, learner_list)))
  
  output <- list(coefs = tmle3_fit$summary, tmle3_fit = tmle3_fit, tmle3_input = tmle3_input)
  class(output) <- c("causalGLMnet", "causalGLM")
  return(output)
}





#' causalRobustGLM
#' Nonparametric robust generalized linear models for causal inference
#' Supports flexible nonparametric conditional average treatment effect (CATE), conditional odds ratio (OR), and conditional relative risk (RR) estimation, 
#' ... where a user-specified working parametric model for the estimand is viewed as an approximation of the true estimand and nonparametrically correct inference is given for these approximations.
#' Specifically, a causal projection of the true estimand onto the working-model is estimated; the parametric model is not assumed correct. 
#' The estimates and inference obtained by `causalRobustGLM` are robust and nonparametrically correct, which comes at a small cost in confidence interval width relative to `causalGLM`.
#' Highly Adaptive Lasso (HAL) (see \code{\link[hal9001]{fit_hal}}), a flexible and adaptive spline regression estimator, is recommended for medium-small to large sample sizes.
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param data A data.frame or matrix containing the numeric values corresponding with the nodes \code{W}, \code{A} and \code{Y}.
#' @param W A character vector of covariates contained in \code{data}
#' @param A A character name for the treatment assignment variable contained in \code{data}
#' @param Y A character name for the outcome variable contained in \code{data} (outcome can be continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "SuperLearner: A stacked ensemble of all of the below that utilizes cross-validation to adaptivelly choose the best learner.
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
#' CATE: Estimate the best parametric approximation of the conditional average treatment effect with \code{\link[tmle3]{Param_npCATE}} assuming it satisfies parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the ATE (Marginal average treatment effect). 
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model.
#' CATT: Estimate the best parametric approximation of the conditional average treatment effect among the treated with \code{\link[tmle3]{Param_npCATE}} assuming it satisfies parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the ATT (Marginal average treatment effect among the treated). 
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model using only the observations with `A=1` (among the treated).
#' OR: Estimate the best parametric approximation of the conditional odds ratio with \code{\link[tmle3]{Param_npOR}} assuming it satisfies parametric model \code{formula}.
#' Specifically, this estimand is the log-likelihood projection of the true conditional odds ratio onto the partially-linear logistic regression model with the true `E[Y|A=0,W]` used as offset.
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
#' @param ... Other arguments to pass to main routine (spCATE, spOR, spRR) 
#' @export 
causalRobustGLM <- function(formula, data, W, A, Y, estimand = c("CATE", "CATT", "OR"),   learning_method = c(  "HAL", "SuperLearner", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"),  cross_fit = FALSE,  sl3_Learner_A = NULL, sl3_Learner_Y = NULL,    weights = NULL,   formula_Y = formula(paste0("~ . + . *", A)),  formula_HAL_Y = paste0("~ . + h(.,", A, ")"), HAL_args_Y = list(smoothness_orders = 1, max_degree = 2, num_knots = c(15,10,1)), HAL_fit_control = list(parallel = F), ... ){
  estimand <- match.arg(estimand)
  learning_method <- match.arg(learning_method)
  data <- as.data.table(data)
  if(!is.null(weights)) {
    data$weights <- weights
  } else {
    data$weights <- 1
  }
  
  if(is.null(sl3_Learner_A)) {
    sl3_Learner_A <- learner_list_A[[learning_method]]
    if(learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
      sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A)
    }
  }
  if(is.null(sl3_Learner_Y)) {
    if(learning_method == "HAL") {
      wrap_in_Lrnr_glm_sp <- FALSE 
      # Allow for formula_HAL
      sl3_Learner_Y <- Lrnr_hal9001$new(formula_HAL = formula_HAL_Y, family = family_list[[estimand]], 
                                        smoothness_orders = HAL_args_Y$smoothness_orders,
                                        max_degree = HAL_args_Y$max_degree,
                                        num_knots = HAL_args_Y$num_knots, fit_control = HAL_fit_control)
      
    } else {
      superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(list(  Lrnr_glmnet$new(formula = formula_Y),  Lrnr_glm$new(formula = formula_Y), Lrnr_gam$new(formula = formula_Y), Lrnr_earth$new(formula = formula_Y) ,
                                                                        Lrnr_ranger$new() ,  Lrnr_xgboost$new(verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =5) ) ), Lrnr_cv_selector$new(loss_function_least_squares))
      learner_list_Y = list(SuperLearner = superlearner_default, glmnet =  Lrnr_glmnet$new(formula = formula_Y), glm = Lrnr_glm$new(formula = formula_Y), gam = Lrnr_gam$new(formula = formula_Y), mars = Lrnr_earth$new(formula = formula_Y) ,
                            ranger = Lrnr_cv$new(Lrnr_ranger$new()), xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,max_depth =5) )), Lrnr_cv_selector$new(loss_squared_error)))
      sl3_Learner_Y <- learner_list_Y[[learning_method]]
    }
    if(learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
      sl3_Learner_Y <- Lrnr_cv$new(sl3_Learner_Y, full_fit = TRUE)
    }
  }
  
  tmle_spec_np <- tmle3_Spec_npCausalGLM$new(formula = formula, estimand = estimand  )
  learner_list <- list(A = sl3_Learner_A, Y = sl3_Learner_Y)
  node_list <- list(weights = "weights", W = W, A = A, Y = Y)
  tmle3_input <- list(tmle_spec_np = tmle_spec_np, data = data, node_list = node_list,  learner_list = learner_list)
  tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec_np, data, node_list, learner_list)))
  
  output <- list(coefs = tmle3_fit$summary,   tmle3_fit = tmle3_fit, tmle3_input = tmle3_input)
  class(output) <- c("causalRobustGLM", "causalGLM")
  return(output)
}
















