


#' npCOXph
#' Nonparametric robust generalized linear models for causal inference and marginal structural models (for some estimands).
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
#' "SuperLearner": A stacked ensemble of all of the below that utilizes cross-validation to adaptivelly choose the best learner.
#' "HAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001} Good for most sample sizes when propertly tuned. See arguments \code{max_degree_Y0W} and \code{num_knots_Y0W}.
#' "glm": Fit nuisances with parametric model. Best for smaller sample sizes (e.g. n =30-100). See arguments \code{glm_formula_A}, \code{glm_formula_T} and \code{glm_formula_T0}.
#' "glmnet": Learn using lasso with glmnet. Best for smaller sample sizes (e.g. n =30-100)
#' "gam": Learn using generalized additive models with mgcv. Good for small-to-medium-small sample sizes.
#' "mars": Multivariate adaptive regression splines with \code{earth}. Good for small-to-medium-small sample sizes.
#' "ranger": Robust random-forests with the package \code{Ranger} Good for medium-to-large sample sizes.
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7. Good for medium-to-large sample sizes.
#' We recommend performing simulations checking 95% CI coverage when choosing learners (especially in smaller sample sizes).
#' Note speed can vary significantly depending on learner choice! 
#' @param estimand Estimand/parameter to estimate. Choices are:
#' `CATE`: Estimate the best parametric approximation of the conditional average treatment effect with \code{\link[tmle3]{Param_npCATE}} using the parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the ATE (Marginal average treatment effect). 
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model.
#' `CATT`: Estimate the best parametric approximation of the conditional average treatment effect among the treated with \code{\link[tmle3]{Param_npCATE}} using the parametric model \code{formula}.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the ATT (Marginal average treatment effect among the treated). 
#' Specifically, this estimand is the least-squares projection of the true CATE onto the parametric working model using only the observations with `A=1` (among the treated).
#' `TSM`: Estimate the best parametric approximation of the conditional treatment-specific mean `E[Y|A=a,W]` for `a` in \code{levels_A}.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the TSM (Marginal treatment-specific mean). 
#' Specifically, this estimand is the least-squares projection of the true TSM onto the parametric working model.
#' `OR`: Estimate the best parametric approximation of the conditional odds ratio with \code{\link[tmle3]{Param_npOR}} using the parametric model \code{formula}.
#' Specifically, this estimand is the log-likelihood projection of the true conditional odds ratio onto the partially-linear logistic regression model with the true `E[Y|A=0,W]` used as offset.
#' `RR`: Projection of the true conditional relative risk onto a exponential working-model using log-linear/poisson regression.
#' Note: if \code{formula} = `~1` then \code{causalRobustGLM} returns a nonparametric and efficient estimator for the marginal relative risk (`E_W[E[Y|A=1,W]]/E_W[E[Y|A=0,W]]``). 
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
#' @param sl3_Learner_T A \code{sl3} Learner object to use to nonparametrically [Y|A,W] with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided. 
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param formula_T Only used if `learning_method %in% c("glm", "earth", "gam", "glmnet")`. A R \code{formula} object that specifies the design matrix to be passed to the Learner specified by learning_method: "glm", "earth", "gam", "glmnet".
#' By default, `formula_T = . + A*.` so that additive learners still model treatment interactions.
#' @param formula_HAL_T A HAL formula string to be passed to \code{\link[hal9001]{fit_hal}}). See the `formula` argument of \code{\link[hal9001]{fit_hal}}) for syntax and example use.
#' @param weights An optional vector of weights to use in procedure.
#' @param HAL_args_T A list of parameters for the semiparametric Highly Adaptive Lasso estimator for E[Y|A,W].
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
npCOXph <- function(formula, data, W, A, Ttilde, Delta, num_bins_t = 20,  learning_method = c(  "HAL", "SuperLearner", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"),  cross_fit = FALSE,  sl3_Learner_A = NULL, sl3_Learner_T = NULL, sl3_Learner_C = NULL,     formula_T = as.formula(paste0("~ . + . *", A)),  formula_HAL_T = paste0("~ . + h(.,", A, ") + h(.,t)"), HAL_args_T = list(smoothness_orders = 1, max_degree = 2, num_knots = c(10,5,1)),  HAL_fit_control = list(parallel = F), delta_epsilon = 0.025, verbose = TRUE, ... ){
   stop("This function is not yet available. ")
  learning_method <- match.arg(learning_method)
  data <- as.data.table(data)
  # if(!is.null(weights)) {
  #   data$weights <- weights
  # } else {
  #   data$weights <- 1
  # }
  nuniq <- length(unique(data[[Ttilde]]))
  if(nuniq > num_bins_t) {
    warning("More than num_bins_t time points. Discretizing...")
    Tt <- data[[Ttilde]]
    bins <-  quantile(ceiling(Tt), seq(0,1, length = num_bins_t))
    Tt <- as.vector(bins[findInterval(Tt, bins, all.inside = TRUE)])
    if(0 %in% bins) {
      Tt <- Tt + 1
    }
    data[[Ttilde]] <- Tt
  }
  
  superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(  Lrnr_glmnet$new(),  Lrnr_glm$new(), Lrnr_gam$new(), Lrnr_earth$new() ,
                                                                         Lrnr_ranger$new() ,  Lrnr_xgboost$new(verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,max_depth =5) ),  full_fit = T) , Lrnr_cv_selector$new(loss_squared_error))
  
  learner_list_A <- list(HAL = Lrnr_hal9001$new(max_degree = 2, smoothness_orders  = 1, num_knots = c(10,3)), SuperLearner = superlearner_default, glmnet =  Lrnr_glmnet$new(), glm = Lrnr_glm$new(), gam = Lrnr_gam$new(), mars = Lrnr_earth$new() ,
                         ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,max_depth =3, eval_metric = "logloss" ), Lrnr_xgboost$new(verbose=0,max_depth =4, eval_metric = "logloss" ), Lrnr_xgboost$new(verbose=0,max_depth =5, eval_metric = "logloss") ), full_fit = TRUE), Lrnr_cv_selector$new(loss_loglik_binomial)))
  
  if(is.null(sl3_Learner_A)) {
    sl3_Learner_A <- learner_list_A[[learning_method]]
    if(learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
      sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A,  full_fit = TRUE)
    }
  }
  binary <- all(Y %in% c(0,1))
  if(is.null(sl3_Learner_T)) {
    if(learning_method == "HAL") {
       
      sl3_Learner_T <- Lrnr_hal9001$new(formula_HAL = formula_HAL_T, family = "binomial", 
                                        smoothness_orders = HAL_args_T$smoothness_orders,
                                        max_degree = HAL_args_T$max_degree,
                                        num_knots = HAL_args_T$num_knots, fit_control = HAL_fit_control)
      
    }  else {
      superlearner_default <- make_learner(Pipeline, Lrnr_cv$new(list(  Lrnr_glmnet$new(formula = formula_T),  Lrnr_glm$new(formula = formula_T), Lrnr_gam$new(), Lrnr_earth$new(formula = formula_T) ,
                                                                        Lrnr_ranger$new() ,  Lrnr_xgboost$new(verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =5) ) , full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error))
      learner_list_Y = list(SuperLearner = superlearner_default, glmnet =  Lrnr_glmnet$new(formula = formula_T), glm = Lrnr_glm$new(formula = formula_T), gam = Lrnr_gam$new(), mars = Lrnr_earth$new(formula = formula_T) ,
                            ranger = Lrnr_cv$new(Lrnr_ranger$new(), full_fit = TRUE), xgboost = make_learner(Pipeline, Lrnr_cv$new( Stack$new(  Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =3 ), Lrnr_xgboost$new(verbose=0,verbose=0,max_depth =4 ), Lrnr_xgboost$new(verbose=0,max_depth =5) ), full_fit = TRUE), Lrnr_cv_selector$new(loss_squared_error)))
      sl3_Learner_T <- learner_list_Y[[learning_method]]
    }
    if(is.null(sl3_Learner_C)) {
      sl3_Learner_C <- sl3_Learner_T
      if(learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
        sl3_Learner_C <- Lrnr_cv$new(sl3_Learner_C, full_fit = TRUE)
      }
    }
    
    if(learning_method %in% c("glm", "glmnet", "mars") && cross_fit) {
      sl3_Learner_T <- Lrnr_cv$new(sl3_Learner_T, full_fit = TRUE)
      
    }
  }
 
  tmle_spec_np <- tmle3_Spec_coxph$new(formula = formula,  delta_epsilon = delta_epsilon, verbose = verbose, treatment_level = 1, control_level = 0)
  learner_list <- list(A = sl3_Learner_A, N = sl3_Learner_T, A_c = sl3_Learner_C )
  node_list <- list(  W = W, A = A, T_tilde = Ttilde, Delta = Delta )
  
  tmle3_input <- list(tmle_spec_np = tmle_spec_np, data = as.data.frame(data), node_list = node_list,  learner_list = learner_list)
  tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec_np, data, node_list, learner_list)))
  
  output <- list(coefs = tmle3_fit$summary,   tmle3_fit = tmle3_fit, tmle3_input = tmle3_input)
  class(output) <- c("npCOXph", "causalGLM")
  return(output)
}




