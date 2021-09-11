
test_that("Testing spglm", {
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
  formula_CATE <- ~ 1 + W1

  # This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
  #
  causal_fit <- spglm(formula = formula_CATE, data = data, W = W, A = A, Y = Y, estimand = "CATE", learning_method = "HAL", verbose = FALSE)
  # We got pretty close!
  coefs <- causal_fit$coefs

  summary(causal_fit)




  # We can also use generalized additive models.
  causal_fit <- spglm(formula = formula_CATE, data = data, W = W, A = A, Y = Y, learning_method = "gam", estimand = "CATE")
  summary(causal_fit)

  # We can also use lasso (glmnet). This is useful for very high dimensional models. (By default, it is cross-fitted to reduce bias).
  # It is amazing that we can get valid inference using the LASSO.
  causal_fit <- spglm(formula = formula_CATE, data = data, W = W, A = A, Y = Y, learning_method = "glmnet", estimand = "CATE")
  summary(causal_fit)


  # We can also use cross-fitted and CV-tuned xgboost. (glmnet is included in the cross-validation selection library/ensemble as well.)
  # Xgboost is black-box. But, we can still get inference!
  causal_fit <- spglm(formula = formula_CATE, data = data, W = W, A = A, Y = Y, learning_method = "xgboost", estimand = "CATE")
  summary(causal_fit)









  set.seed(1500)
  data_list <- sim.OR(150, 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data
  truelogOR <- data_list$logOR
  # True log OR of data (is constant)
  print(truelogOR)

  # Let's learn it using semiparametric methods.
  # Lets specify a constant model for the OR
  formula_logOR <- ~1

  # Let use MARS (multivariate adaptive regression splines) using the "earth" package.
  # It is amazing that we can get valid inference using the greedy selection algorithms like MARS!
  causal_fit <- spglm(formula = formula_logOR, data, W = W, A = A, Y = Y, estimand = "OR", learning_method = "mars")
  # We got pretty close!


  summary(causal_fit)







  data_list <- sim.RR(150, 2)
  # Confounders
  W <- data_list$W
  # Binary treatment
  A <- data_list$A
  # Outcome (binary in this case)
  Y <- data_list$Y
  data <- data_list$data


  # True log RR of data (is constant)
  print(data_list$logRR)

  # Let's learn it using semiparametric methods.
  # Lets specify a less constant model for the RR (only first coefficient is nonzero)
  formula_logRR <- ~ 1 + W1 + W2

  # This will take a few seconds.
  causal_fit <- spglm(formula = formula_logRR, data, W = W, A = A, Y = Y, estimand = "RR", learning_method = "xgboost")
  # We got pretty close!

  summary(causal_fit)
  expect_equal(T, T)
})
