
 
 
# causalglm : interpretable and model-robust causal inference for heterogeneous treatment effects


In the search for answers to causal questions, assuming parametric models can be dangerous. With a seemingly small amount of confounding and model misspecificaton, they can give biased answers. One way of mitigating this challenge is to only parametrically model the feature of the data-generating distribution that you care about. It is not even necessary to assume your parametric model is correct! Instead, view it as a "working" or "approximation" model and define your estimand as the best causal approximation of the true nonparametric estimand with respect to your parametric working model. This allows for causal estimates and robust inference under no parametric assumptions on the functional form of any feature of the data generating distribution. Alternatively, you can assume a semiparametric model that only assumes the parametric form of the relevant part of the data distribution is correct. Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question.  It is in fact possible to get robust and efficient inference for causal quantities using machine-learning. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.


This package fully utilizes the powerful `tlverse/tmle3` generalized targeted learning framework as well as the machine-learning frameworks `tlverse/sl3` and `tlverse/hal9001`. We recommend taking a look at these packages and the rest of the `tlverse`! 

For in progress theoretical details and methods descriptions, see the writeup `causalglm.pdf` in the "writeup" folder.

For a walk-through guide with example code, check out `vignette.Rmd` in the "vignette" folder.

## Installing
NOTE: This package is actively in development and is subject to continuous change. Feel free to contact me personally or through the Issues tab.

 
 
To install this package, install the devtools CRAN package and run:

``` r
if(!require(devtools)) {
  install.packages(devtools)
}
devtools::install_github("tlverse/causalglm")
```

This package also requires the following github packages. Make sure to update the version of these packages upon installation.
``` r
devtools::install_github("tlverse/hal9001@master")
devtools::install_github("tlverse/tmle3@general_submodels_devel")
devtools::install_github("tlverse/sl3@Larsvanderlaan-formula_fix")
```


## What is causalglm?


`causalglm` is a package for interpretable and robust causal inference for heterogeneous treatment effects using generalized linear working model and machine-learning. Unlike parametric methods based on generalized linear models and semiparametric methods based on partially-linear models, the methods implemented in `causalglm` do not assume that any user-specified parametric models are correctly specified. That is, `causalglm` does not assume that the true data-generating distribution satisfies the parametric model. Instead, the user-specified parametric model is viewed as an approximation or "working model", and an interpretable estimand is defined as a projection of the true conditional treatment effect estimand onto the working model. Moreover, `causalglm`, unlike `glm`, only requires a user-specified parametric working model for the causal estimand of interest. All nuisance components of the data-generating distribution are left unspecified and data-adaptively learned using machine-learning. Thus, `causalglm` provides not only nonparametrically robust inference but also provides estimates (different from `glm`) for causally interpretable estimands that maximally adjust for confounding. To allow for valid inference with the use of variable-selection and machine-learning, Targeted Maximum Likelihood Estimation (van der Laan, Rose, 2011) is employed.  

The statistical data-structure used throughout this package is `O = (W,A,Y)` where `W` represents a random vector of baseline (pretreatment) covariates/confounders, `A` is a binary, categorical or continuous treatment assignment, and `Y` is some outcome variable. For marginal structural models, we also consider a subvector `V` of `W` that represents a subset of baseline variables that are of interest.

### Estimands supported by `causalglm`

`causalglm` supports causal working-model-based estimands for the

1. Conditional average treatment effect (CATE) for arbitrary outcomes: `E[Y|A=a,W] - E[Y|A=0,W]` (categorical and continuous treatments)
2. Conditional odds ratio (OR) for binary outcomes: `{P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)}` (binary treatments and continuous treatments)
3. Conditional relative risk (RR) for binary, count or nonnegative outcomes: `E[Y|A=a,W]/E[Y|A=0,W]` (categorical and continuous treatments)
4. Conditional treatment-specific mean (TSM) : `E[Y|A=a,W]` (categorical treatments)
5. Conditional average treatment effect among the treated (CATT) : the best approximation of `E[Y|A=a,W] - E[Y|A=0,W]` based on a user-specified formula/parametric model among the treated (i.e. observations with `A=a`) (categorical treatments)

All methods support binary treatments. Most methods support categorical treatments. And, continuous treatments are only supported for the `CATE`, `OR` and `RR` through `contglm`. Each method allows for arbitrary user-specified parametric working models for the estimands. For binary and categorical treatments, the working model used for all estimands is of the form `E[Y|A=a,W] - E[Y|A=0,W] = formula(W)` where `formula` is specified by the user. For continuous treatments (only supported by `contglm`), the working models used are `CATE(a,W) := E[Y|A=a,W] - E[Y|A=0,W] = 1(a > 0) * formula_binary(W) + a * formula_continuous(W)`, `log OR(a,W) := log {P(Y=1|A=a,W)/P(Y=0|A=a,W) } - log {P(Y=1|A=0,W)/P(Y=0|A=0,W) }= 1(a > 0) * formula_binary(W) + a * formula_continuous(W)`, and `log RR(a,W) := log E[Y|A=a,W]  - log E[Y|A=0,W] = 1(a > 0) * formula_binary(W) + a * formula_continuous(W)` Since estimates and inference are provided for the best approximation of the estimand by these working models, the outputs of `contglm` can be viewed as the best linear approximation to the continuous treatment estimand.

causalglm also supports the following working marginal structural model estimands:
 
1. Working marginal structural models for the CATE: `E[CATE(W)|V] := E[E[Y|A=a,W] - E[Y|A=0,W]|V]` (categorical treatments)
2. Working marginal structural models for the RR: `E[E[Y|A=a,W]|V]/E[E[Y|A=0,W]|V]` (categorical treatments)
3. Working marginal structural models for the TSM : `E[E[Y|A=a,W]|V]` (categorical treatments)
4. Working marginal structural models for the CATT : `E[CATE(W)|V, A=a] := E[E[Y|A=a,W] - E[Y|A=0,W]|V, A=a]` (categorical treatments)
 

### Methods provided by `causalglm`
causalglm consists of 2 main functions: 
 
1. `npglm` for robust nonparametric estimation and inference for user-specified working models for the `CATE`, `CATT`, `TSM`, `RR` or `OR`
2. `contglm` for robust nonparametric estimation and inference for user-specified working models for the `CATE`, `OR` and `RR` as a function of a continuous or ordered numeric treatment.

And 3 more specialized functions:

4. `msmglm` for robust nonparametric estimation and inference for user-specified working marginal structural models for the `CATE`, `CATT`, `TSM` or `RR`
5. `spglm` for semiparametric estimation and inference for correctly specified parametric models for the `CATE`, `RR` and `OR`
6. `causalglmnet` for semiparametric estimation and inference with high dimensional confounders `W` (a custom wrapper function for spglm focused on big data where standard ML may struggle)

For most user applications with discrete treatments, `npglm` suffices. For continuous treatments, users may use `contglm`.
 

 
### Outputs provided by all `causalglm` methods

The outputs of the methods include:

1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients
4. Individual-level treatment-effect predictions and 95\% confidence (prediction) intervals can be extracted with the `predict` function and argument `data`.
5. Plotting with `plot_msm` for objects returned by `msmglm`.


### Which method to use?
A rule of thumb for choosing between these methods is as follows:

1. Use `npglm` if you believe your parametric model for the treatment effect estimand is a good approximation
2. Use `msmglm` if you want to know how the treatment effect is causally affected by one or a number of variables `V` (fully adjusting for the remaining variables `W`) (or to learn univariate confounder-adjusted variable importance measures!)
3. Use `contglm` if your treatment is continuous or ordered and you are interested in the treatment effect per unit dose.
4. Use `causalglmnet` if the variables `W` for which to adjust are (very) high dimensional.
5. Use `spglm` if you believe your parametric model for the treatment effect estimand is correct (not recommended)
 

`msmglm` deals with marginal structural models for the conditional treatment effect estimands. This method is useful if you are only interested in modeling the causal treatment effect as a function of a subset of variables `V` adjusting for all the available confounders `W` that remain. This allows for parsimonious causal modeling, still maximally adjusting for confounding. This function can be used to understand the causal variable importance of individual variables (by having `V` be a single variable) and allows for nice plots (see `plot_msm`). `contglm` is a version of `npglm` that provides inference for working-model-based estimands for conditional treatment effects of continuous or ordered treatments.  
 


 
### User-friendly interface
`causalglm` has a minimalistic yet still quite flexible front-end. Check out the vignette to see how to use it! The necessary arguments are: 
1. `formula`: a formula object for the `CATE`, `OR`, or `RR` (also `TSM`, `CATT` for `npglm`)
2. `data`: a data.frame containing the data
3. `W`, `A`, `Y`:  character vectors that store the variable names for the baseline variables, treatment variable and outcome variable.
4. `estimand`: the choice of estimand, `"CATE"`, `"OR"`, `"RR"` (`"TSM"`, `"CATT"` for `npglm`)
 
Other useful arguments are:

5. `learning_method`: The machine-learning algorithm used to learn the nuisance function (default is HAL. See vignette and documentation)
6. `treatment_level` for specifying what value of a possibly categorical `A` is the treatment of interest (`A=a`).
7. `control_level` for specifying what value of a possibly categorical `A` is the control of interest (`A=0`).
8. `formula_Y` for specifying the design matrix passed to glm, glmnet, earth or gam.
9. `formula_HAL_Y` for specifying a custom functional form for the Highly Adaptive Lasso estimator (HAL) 
 

### Conditional average treatment effect estimation


``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A * (1 + W + 2*W^2) + sin(5*W), sd = 0.3)
data <- data.frame(W,A,Y)

 

# Nonparametric
formula <- ~ poly(W, degree = 2, raw = FALSE)
output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "CATE",
    learning_method = "HAL",
    formula_HAL_Y = ~ h(W) + h(W,A), # Optional
    verbose = FALSE
  )

summary(output) 
head(predict(output, data = data))
 
 # Inference for best linear approximation
formula <- ~ 1 + W
output <-
  npglm(
    formula,
    output, # Reuse ML fits of previous output by passing in the output as the data argument
    estimand = "CATE",
    verbose = FALSE
  )

summary(output) 
head(predict(output, data = data))


# Learn a marginal structural model for the CATE as a function of `V`.
formula <- ~ poly(W, degree = 2, raw = FALSE)
output <-
  msmglm(
    formula,
    data,
    V = "W",
    W = "W", A = "A", Y = "Y",
    estimand = "CATE",
    learning_method = "glmnet",
    formula_Y = ~ . ^2, # Optional formula for outcome model used to for glmnet estimation
    verbose = FALSE
  )

summary(output) 
plot_msm(output)

 
```

### CATE with categorical treatments
(Same idea for CATT, TSM and RR)
``` r
n <- 500
W <- runif(n, min = -1, max = 1)
A <- rbinom(n, size = 1, prob = 0.66 * plogis(W))
A[A == 1] <- 2
A[A == 0] <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A * (1 + W) + W , sd = 0.5)
data <- data.table(W, A, Y)

# Model is E[Y|A=treatment_level, W] -  E[Y|A=control_level, W] = formula(W)

output_init <- npglm(
  formula = ~ 1 + W,
  data,
  W = "W", A = "A", Y = "Y",
  estimand = "CATE",
  learning_method = "mars",
  treatment_level = 1,
  control_level = 0
)

summary(output_init)


# Reuse fits
output <- npglm(
  ~ 1 + W,
  data = output_init ,
  estimand = "CATE",
  treatment_level = 2,
  control_level = 0
)


summary(output)

```


### CATE with continuous treatments

Change the `estimand` argument to estimate the RR and OR.
``` r
n <- 500
W <- runif(n, min = -1, max = 1)
Abinary <- rbinom(n , size = 1, plogis(W))
A <- rgamma(n, shape = 1, rate = exp(W))
A <- A * Abinary
Y <- rnorm(n, mean =   (A > 0) + A * (1 + W) + W , sd = 0.5)
data <- data.table(W, A, Y)

# Model is CATE(A,W) = formula_binary(W) 1(A > 0) + A * formula_continuous(W)

out <- contglm(
  formula_continuous = ~ 1 + W,
  formula_binary = ~ 1,
  data = data,
  W = "W", A = "A", Y = "Y",
  estimand = "CATE",
  learning_method = "gam",
  formula_Y = ~ . + . * A # Optional: Feed design matrix based on formula to gam
)

summary(out)

# The CATE predictions are now a function of `A`
head(predict(out))

```


### Conditional odds ratio estimation
Note a log-linear working model is used for the conditional odds ratio.
As a consequence, the parametric working model specified by `formula` is actually a working model for the log conditional odds ratio.

``` r
# odds ratio
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n, size =  1, prob = plogis(A + A * W + W + sin(5 * W)))
data <- data.frame(W, A, Y)
# Nonparametric robust inference
output <-
  npglm(
    ~1+W,
    data,
    W = c("W"), A = "A", Y = "Y",
    estimand = "OR",
    learning_method = "xgboost" # Default xgboost ensemble
  )
summary(output)

```

 


### Conditional relative risk/treatment-effect estimation
 Note a log-linear working model is used for the conditional relative risk.
As a consequence, the parametric working model specified by `formula` is actually a working model for the log conditional relative risk.

 

``` r
# relative risk
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp( A * (1 + W + 2*W^2)  + sin(5 * W)))
data <- data.frame(W, A, Y)
formula = ~ poly(W, degree = 2, raw = TRUE) 
# Nonparametric robust inference
output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "RR",
    estimand = "OR",
    learning_method = "SuperLearner", # Default SuperLearner ensemble using all methods
    formula_Y = ~.^2, # formula used for glm, glmnet, gam and earth in SuperLearner
    verbose = FALSE
  )
summary(output)

# Semiparametric inference
output <-
  spglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "RR",
    verbose = FALSE
  )
summary(output)

# Learn a working marginal structural model with cool plots
output <-
  msmglm(
    formula,
    data,
    V = "W",
    W = "W", A = "A", Y = "Y",
    estimand = "RR",
    learning_method = "HAL",
    formula_HAL_Y = ~ h(.) + h(.,.),
    verbose = FALSE
  )
summary(output)
 
```

### More conditional treatment effect estimation with npglm and msmglm (the CATT and TSM)

``` r
n <- 250
W1 <- runif(n, min = -1, max = 1)
W2 <- runif(n, min = -1, max = 1)
A <- rbinom(n, size = 1, prob = plogis((W1 + W2  )/3))
Y <- rnorm(n, mean = A * (1 + W1 + 2*W1^2) + sin(4 * W2) + sin(4 * W1), sd = 0.3)
data <- data.frame(W1, W2,A,Y)
# CATE 
formula = ~ poly(W1, degree = 2, raw = TRUE)
output <- npglm(formula,
      data,
      W = c("W1", "W2"), A = "A", Y = "Y",
      estimand = "CATE")
summary(output)
# CATT, lets reuse the fit (see vignette)
output <- npglm(formula,
      output,
      estimand = "CATT")
summary(output)
# TSM, note this provides a list of npglm objects for each level of `A`.
outputs <- npglm(formula,
      output,
      estimand = "TSM")
summary(outputs[[1]])
summary(outputs[[2]])


formula = ~ poly(W1, degree = 2, raw = TRUE)
output <- msmglm(formula,
      data,
      V = "W1",
      W = c("W1", "W2"), A = "A", Y = "Y",
      estimand = "CATE")
summary(output)
plot_msm(output)
```

## Using tlverse/sl3 to specify custom learners.
Custom learners can be specified for the conditional means of `A` and `Y` using the `sl3_Learner_A` and `sl3_Learner_Y` arguments. 
``` r
n <- 250
W1 <- runif(n, min = -1, max = 1)
W2 <- runif(n, min = -1, max = 1)
A <- rbinom(n, size = 1, prob = plogis((W1 + W2  )/3))
Y <- rnorm(n, mean = A * (1 + W1 + 2*W1^2) + sin(4 * W2) + sin(4 * W1), sd = 0.3)
data <- data.frame(W1, W2,A,Y)

# Use sl3
library(sl3)
# All sorts of options
lrnr_A <- Lrnr_glm$new(formula = ~ .^2)
lrnr_Y <- Lrnr_glmnet$new(formula = ~ .^2)
lrnr_Y <- Lrnr_xgboost$new(max_depth  = 5)
lrnr_Y <- Lrnr_gam$new()
lrnr_Y <- Lrnr_earth$new()
lrnr_stack <- make_learner(Stack, Lrnr_glmnet$new(), Lrnr_glm$new(), Lrnr_gam$new(), Lrnr_xgboost$new()) # Add as many learners as you want in the stack
lrnr_sl <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new()) # CV selection

# CATE estimation
formula = ~ poly(W1, degree = 2, raw = TRUE)
output <- npglm(formula,
      data,
      W = c("W1", "W2"), A = "A", Y = "Y",
      estimand = "CATE",
      sl3_Learner_A = lrnr_A,
      sl3_Learner_Y = lrnr_Y)
summary(output)

```

## Learn working marginal structural models for conditional treatment effects with `msmglm`

The `msmglm` function is a wrapper for `npglm` that focuses purely on working marginal structural model estimation and it has some convenient plotting features. `V` can be multivariate but plotting is only supports for the univariate case. `estimand = "CATE"` corresponds with the estimand `E[CATE(W)|V]`, `estimand = "CATT"` corresponds with the estimand `E[CATE(W)|V, A=1]`, `estimand = "TSM"` corresponds with the estimand `E[E[Y|A=a,W]|V]`, and  `estimand = "RR"` corresponds with the estimand `E[E[Y|A=1,W]|V] / E[E[Y|A=0,W]|V]`.  The intercept model reduces to nonparametric efficient estimators for the marginal ATE, ATT, TSM and RR respectively. Note that formula should only depend on `V`. This method is useful to study confounder-adjusted associations between a continuous variable (i.e. `V`) and the treatment effect or outcome `Y`. 


``` r

n <- 250
V <- runif(n, min = -1, max = 1)
W <- runif(n, min = -1, max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))

# CATE
Y <- rnorm(n, mean = A * (1 + V + 2*V^2) + W + V + sin(5 * W), sd = 0.5)
data <- data.frame(V,W, A, Y)
formula_msm <- ~ poly(V, degree = 2, raw = TRUE) # A second degree polynomial
output <-
  msmglm(
    formula_msm,
    data,
   V = "V",
    W = c("V","W"), A = "A", Y = "Y",
    estimand = "CATE",
    learning_method = "glm",
    formula_Y = ~ . + . * A,
    verbose = FALSE
  )

summary(output)
plot_msm(output)

# CATT
Y <- rnorm(n, mean = A * (1 + V + 2*V^2) + W + V + sin(5 * W), sd = 0.5)
data <- data.frame(V,W, A, Y)
formula_msm <- ~ poly(V, degree = 2, raw = TRUE) 
output <-
  msmglm(
    formula_msm,
    output,
    V = "V",
    estimand = "CATT",
    verbose = FALSE
  )

summary(output)
plot_msm(output)

# TSM
Y <- rnorm(n, mean = A * (1 + V + 2*V^2) + W + V, sd = 0.5)
data <- data.frame(V,W, A, Y)
formula_msm <- ~ poly(V, degree = 2, raw = TRUE) 
output <-
  msmglm(
    formula_msm,
    output,
    V = "V",
    estimand = "TSM",
    verbose = FALSE
  )

summary(output[[1]])
summary(output[[2]])
 plot_msm(output[[1]])

# RR
Y <- rpois(n, lambda = exp( A * (1 + V + 2*V^2)  + sin(5 * W)))
data <- data.frame(V,W, A, Y)
formula_msm <- ~ poly(V, degree = 2, raw = TRUE) 
output <-
  msmglm(
    formula_msm,
    data,
    V = "V",
    W = c("V","W"), A = "A", Y = "Y",
    estimand = "RR",
    verbose = FALSE
  )

summary(output)
plot_msm(output)
```

## Semiparametric inference for high dimensional generalized linear models with causalglmnet (the LASSO): CATE, OR, and RR
For high dimensional W, you can use the wrapper function `causalglmnet` which runs `spglm` using a custom glmnet-LASSO learner for estimation. This allows for robust and fast estimation in high dimensional settings where conventional machine-learning algorithms may struggle. Cross-fitting can be performed to reduce bias. This method can be viewed as an adaptive version of "glm" in that confounders/variables to adjust for are adaptively selected using the LASSO, while still allowing for asymptotically correct post-selection inference. 

``` r
n <- 200
W <- replicate(100,runif(n, min = -1,  max = 1))
colnames(W) <- paste0("W", 1:100)
beta <- runif(10, -1, 1)/20
A <- rbinom(n, size = 1, prob = plogis(W[,10*(1:10)] %*% beta))

# CATE
Y <- rnorm(n, mean = A + W[,10*(1:10)] %*% beta, sd = 0.5)
data <- data.frame(W, A, Y)

formula = ~ 1  
output <-
  causalglmnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "CATE" ,
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))

# OR
Y <- rbinom(n, size =  1, prob = plogis( A + W[,10*(1:10)] %*% beta))
data <- data.frame(W, A, Y)

formula  = ~ 1  
output <-
  causalglmnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "OR" ,
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))

# RR
Y <- rpois(n, lambda = exp( A + W[,10*(1:10)] %*% beta))
data <- data.frame(W, A, Y)

formula  = ~ 1  
output <-
  causalglmnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "RR" ,
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))
```

## Need a new or specialized method? Questions? Suggestions?

Any confusion? Questions? Don't know which method to use? None of the methods handle your problem? Need a custom/specialized method?
  
  Just send me a message. I would enjoy connecting with you and be happy to develop and implement a new method to add to this package. Please contact me at vanderlaanlars@yahoo.com.

 
## References (see writeup)

Here is a summary of the methods:

| ----method (Top) / feature (Left) ---| causalglm (as a whole) |npglm | msmglm | spglm | contglm | causalglmnet | glm |
|--------------------------------------|-----------------------|------|---------|-------|---------|--------------|-----|
| Semiparametric inference       |   Y  |      Y  |    Y    |   Y   |    Y    |     Y        |  N  |
| Robust nonparametric working-model-based inference |   Y  |      Y  |    Y    |   N   |    Y    |     N        |  N  |            
| Binary treatment               |   Y  |      Y  |    Y    |   Y   |    N    |     Y        |  Y  |
| Categorical treatment          |   Y  |       Y  |    Y    |   N   |    N    |     N        |  Y  |
| Continuous or ordered treatment|   Y  |       N  |    N    |   N   |    Y    |     N        |  Y  |
| Marginal structural working models     |   Y  |       Y  |    Y    |   N   |    N    |     N        |  N  |
| Interpretable estimates       |   Y  |      Y  |    Y    |   Y   |    Y    |     Y        |  Y  |
| Causal (unconfounded) estimates  under model mispecification |  Y  |   Y  |    Y    |   N   |    Y    |     N        |  N  |
| Supports inference with machine-learning and variable selection |   Y  |   Y  |    Y    |   Y   |    Y    |     Y        |  N  |
| Inference with High dimensional confounders  |   Y  |    Y  |    Y    |   Y   |    Y    |     Y*1        |  N  |
| CATE |   Y   | Y | Y | Y |Y |Y | Y|
| OR |   Y  | Y | Y | Y |Y |Y | Y|
| RR |   Y   | Y | Y | Y |Y |Y | N*2 |
| TSM |   Y   | Y | Y | N |N |N | Y |
| CATT |   Y   | Y | Y | N |N |N |Y|
| p-values and confidence intervals |   Y  |  Y | Y | Y |Y |Y | Y|
| Individual treatment effects with confidence intervals |   Y  |  Y | Y | Y |Y |Y | N|

*1: All methods but glm support the LASSO for estimation of all nuisance parameters and can thus be used in very high dimensions. 
However, causalglmnet uses a customized LASSO learner that should perform better than the other methods in high dimensions.


*2: glm only supports correct inference for the RR when outcome error distribution is poisson. causalglm makes no assumptions on the error distribution and works for arbitrary binary, count and nonnegative outcomes



 
 

