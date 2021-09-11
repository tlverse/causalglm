
context("Testing causalglmnet")

set.seed(1500)
data_list <- sim.CATE(n=100, p=50)
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
formula_CATE <- ~ 1

# This will take a few seconds learning_method = HAL is a regularized sparse/smoothness adaptive regression spline algorithm. To make it faster change: learning_method or max_degree_Y0W or num_knots_Y0W or parallel.
# 
doMC::registerDoMC(5)
causal_fit <- causalglmnet(formula = formula_CATE, data = data, W  = W, A = A, Y = Y,  estimand = "CATE", cross_fit = TRUE, constant_variance_CATE = TRUE, parallel = TRUE)
# We got pretty close!
coefs <- causal_fit$coefs

# The intercept model is actually a nonparametric estimate of the ATE
summary(causal_fit)




set.seed(1500)
data_list <- sim.OR(n=100, p=50)
# Confounders
W <- data_list$W
# Binary treatment
A <- data_list$A
# Outcome (binary in this case)
Y <- data_list$Y
data <- data_list$data
trueOR <- data_list$beta_logRR
# True conditional OR of data (is constant)
 

# Let's learn it using semiparametric methods.
# Lets specify a constant model for the OR
formula_OR <- ~ 1


doMC::registerDoMC(5)
causal_fit <- causalglmnet(formula = formula_CATE, data = data, W  = W, A = A, Y = Y,  estimand = "OR", cross_fit = TRUE, constant_variance_CATE = TRUE, parallel = TRUE)
# We got pretty close!
coefs <- causal_fit$coefs

# The intercept model can be interpreted as a population-average of the conditional odds ratio.
summary(causal_fit)


set.seed(1500)

data_list <- sim.RR(100, 50)
# Confounders
W <- data_list$W
# Binary treatment
A <- data_list$A
# Outcome (binary in this case)
Y <- data_list$Y
data <- data_list$data


# True log RR of data (is constant)
 

# Let's learn it using semiparametric methods.
# Lets specify a less constant model for the RR (only first coefficient is nonzero)
formula_logRR <- ~ 1 + W1 + W2 
doMC::registerDoMC(5)
# This will take a few seconds. 
causal_fit <- causalglmnet(formula = formula_logRR, data, W  = W, A = A, Y = Y,  estimand = "RR", parallel = TRUE, verbose = TRUE   )
# We got pretty close!

summary(causal_fit) 