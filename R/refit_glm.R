
#' internal use to reuse fits
refit_glm <- function(fit_object, formula, estimand = fit_object$estimand, treatment_level = NULL, control_level = NULL, verbose = TRUE) {
  likelihood <- fit_object$tmle3_fit$likelihood$initial_likelihood
  tmle3_input <- fit_object$tmle3_input
  data <- tmle3_input$data
  node_list <- tmle3_input$node_list
  learner_list <- NULL
  delta_epsilon <- tmle3_input$delta_epsilon

  if (inherits(fit_object, "npglm") || inherits(fit_object, "msmglm")) {
    tmle_spec <- tmle3_Spec_npCausalGLM$new(likelihood_override = likelihood, formula = formula, estimand = estimand, delta_epsilon = delta_epsilon, verbose = verbose, treatment_level = treatment_level, control_level = control_level)
  } else {
    stop("Reusing fits is not supported for spglm")
  }
  tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec, data, node_list, learner_list)))

  return(tmle3_fit)
}
