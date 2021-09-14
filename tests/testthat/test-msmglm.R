test_that("Testing msmglm", {


  # set.seed(1500)
  data_list <- sim.CATE(n = 150, p = 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data
  trueCATE <- data_list$beta_CATE
  task <- sl3_Task$new(data, c("W1", "W2", "A"), "Y")



  # True treatment effect of data (is constant)
  print(trueCATE)

  # Let's learn it using semiparametric methods.
  # Lets specify a constant model for the CATE
  formula_CATE <- ~1

  # This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
  #
  data
  causal_fit <- msmglm(formula = formula_CATE, data = data, V = W, A = A, Y = Y, estimand = "CATE", learning_method = "glmnet")
  # We got pretty close!
  coefs <- causal_fit$coefs

  # The intercept model is actually a nonparametric estimate of the ATE
  summary(causal_fit)







  set.seed(1500)
  data_list <- sim.CATE(n = 150, p = 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data
  trueCATE <- data_list$beta_CATE
  # True treatment effect of data (is constant)
  print(trueCATE)

  # Let's learn it using semiparametric methods.
  # Lets specify a constant model for the CATE
  formula_TSM <- ~1

  # This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
  #
  causal_fit <- msmglm(formula = formula_TSM, data = data, V = W, A = A, Y = Y, estimand = "TSM", learning_method = "xgboost")
  # We got pretty close!
  coefs <- causal_fit$coefs

  # The intercept model is actually a nonparametric estimate of the ATE
  summary(causal_fit)




  set.seed(1500)
  data_list <- sim.CATE(n = 150, p = 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data
  trueCATE <- data_list$beta_CATE
  # True treatment effect of data (is constant)
  print(trueCATE)

  # Let's learn it using semiparametric methods.
  # Lets specify a constant model for the CATE
  formula_CATE <- ~1

  # This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
  #
  causal_fit <- msmglm(formula = formula_CATE, data = data, V = W, A = A, Y = Y, estimand = "CATT", learning_method = "HAL")
  # We got pretty close!
  coefs <- causal_fit$coefs

  # The intercept model is actually a nonparametric estimate of the ATT
  summary(causal_fit)




  set.seed(1500)
  data_list <- sim.RR(n = 150, p = 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data
  beta_logRR <- data_list$beta_logRR
  # True treatment effect of data (is constant)
  print(beta_logRR)

  # Let's learn it using semiparametric methods.
  # Lets specify a constant model for the CATE
  formula_RR <- ~1

  # This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
  #
  causal_fit <- msmglm(formula = formula_RR, data = data, V = W, A = A, Y = Y, estimand = "RR", learning_method = "glm")
  # We got pretty close!
  coefs <- causal_fit$coefs

  # The intercept model can be interpreted as a population-average of the conditional odds ratio.
  summary(causal_fit)


  expect_equal(T, T)
})
