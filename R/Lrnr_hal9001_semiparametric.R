#' Scalable Highly Adaptive Lasso (HAL) semiparametric
#' @importFrom R6 R6Class
#' @importFrom origami folds2foldvec
#' @importFrom stats predict quasibinomial
#'
 
Lrnr_hal9001_semiparametric <- R6Class(
  classname = "Lrnr_hal9001_semiparametric",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(formula_sp, family = NULL, interaction_variable = "A", ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),
    .train = function(task) {
      args <- self$params


      formula_sp <- args$formula_sp
      trt <- args$interaction_variable
      A <- task$data[[trt]]
      W <- as.matrix(task$X)
      W <- W[, setdiff(task$nodes$covariates, trt), drop = F]

      V <- model.matrix(formula_sp, as.data.frame(W))

      outcome_type <- self$get_outcome_type(task)
      args$X <- as.matrix(W)
      args$X_unpenalized <- as.matrix(A * V)
      args$Y <- as.numeric(as.character(outcome_type$format(task$Y)))

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family()
      }

      if (!any(grepl("fit_control", names(args)))) {
        args$fit_control <- list()
      }
      args$fit_control$foldid <- origami::folds2foldvec(task$folds)

      if (task$has_node("id")) {
        args$id <- task$id
      }

      if (task$has_node("weights")) {
        args$fit_control$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      # fit HAL, allowing glmnet-fitting arguments
      other_valid <- c(
        names(formals(glmnet::cv.glmnet)), names(formals(glmnet::glmnet))
      )

      fit_object <- sl3:::call_with_args(
        hal9001::fit_hal, args,
        other_valid = other_valid
      )

      return(fit_object)
    },

    .predict = function(task = NULL) {
      args <- self$params
      formula_sp <- args$formula_sp
      trt <- args$interaction_variable
      A <- task$data[[trt]]
      W <- as.matrix(task$X)
      W <- W[, setdiff(task$nodes$covariates, trt)]
      V <- model.matrix(formula_sp, as.data.frame(W))
      X_unpenalized <- as.matrix(A * V)
      predictions <- stats::predict(
        self$fit_object,
        new_data = as.matrix(W),
        new_X_unpenalized = X_unpenalized
      )
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },

    .required_packages = c("hal9001", "glmnet")
  )
)
