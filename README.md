# causalglm : interpretable and robust causal inference for heterogeneous treatment effects


It is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for answers to causal questions, assuming parametric models can be dangerous. With even a seemingly small amount of confounding and misspecificaton, they can give biased answers. One way of mitigating this challenge is to only parametrically model the feature of the data-generating distribution that you care about. That is, use a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.


This package fully utilizes the powerful `tlverse/tmle3` generalized targeted learning framework as well as the machine-learning frameworks `tlverse/sl3` and `tlverse/hal9001`. We recommend taking a look at these packages and the rest of the `tlverse`! 

For theoretical details and methods descriptions, see the writeup `causalglm.pdf` in the "writeup" folder.

For a walk-through guide with example code, check out `vignette.Rmd` in the "vignette" folder.

## Installing
NOTE: This package is actively in development and is subject to continuous change. If you are unable to install it due to errors, try again in a day or two. Also, feel free to contact me personally or through the Issues tab.

 
 


To install this package, install the devtools CRAN package and run:

``` r
if(!require(devtools)) {
  install.packages(devtools)
}
devtools::install_github("tlverse/causalglm")
```

This package also requires the following github packages which should be installed automatically.
``` r
devtools::install_github("tlverse/hal9001@devel")
devtools::install_github("tlverse/tmle3@general_submodels_devel")
devtools::install_github("tlverse/sl3@Larsvanderlaan-formula_fix")
```


## What is causalglm?


causalglm is an R package for robust generalized linear models and interpretable causal inference for heterogeneous (or conditional) treatment effects. Specifically, causalglm very significantly relaxes the assumptions needed for useful causal estimates and correct inference by employing semi and nonparametric models and adaptive machine-learning through targeted maximum likelihood estimation (TMLE) (van der Laan, Rose, 2011), while still robustly utilizing user-specified parametric forms for the conditional estimands and thus allowing for interpretable inference. Because of this, `causalglm` methods can (and often do in real-world settings) significantly reduce bias relative to conventional fully parametric methods like `glm`, and the cost in variance/confidence-interval-width is negligible. As another consequence, causalglm can be used in high dimensional settings and provides valid inference even when adaptive variable selection is used with, for example, the LASSO. See the writeup causalglm.pdf for a more theoretical overview of the methods implemented in this package. Also see the vignette for an overview of all the functionalities of causalglm. 

Throughout the development of causalglm, we greatly focused on making this package and all its methods as easy to use and understand as possible. Our goal is that this package can be used by a wide audience including amateurs, practioners, statisticians, causal-inference experts and nonexperts.

The statistical data-structure used throughout this package is `O = (W,A,Y)` where `W` represents a random vector of baseline (pretreatment) covariates/confounders, `A` is a binary, categorical or continuous treatment assignment, and `Y` is some outcome variable. For marginal structural models, we also consider a subvector `V \subset W` that represents a subset of baseline variables that are of interest.

The estimands supported by causalglm are

1. Conditional average treatment effect (CATE) for arbitrary outcomes: `E[Y|A=a,W] - E[Y|A=0,W]` (categorical and continuous treatments)
2. Conditional odds ratio (OR) for binary outcomes: `{P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)}` (binary treatments)
3. Conditional relative risk (RR) for binary, count or nonnegative outcomes: `E[Y|A=a,W]/E[Y|A=0,W]` (categorical treatments)
4. Conditional treatment-specific mean (TSM) : `E[Y|A=a,W]` (categorical treatments)
5. Conditional average treatment effect among the treated (CATT) : the best approximation of `E[Y|A=a,W] - E[Y|A=0,W]` based on a user-specified formula/parametric model among the treated (i.e. observations with `A=a`) (categorical treatments)

All methods support binary treatments. Most methods support categorical treatments. And, continuous treatments are only supported for the `CATE` through `contglm`. Each method allows for arbitrary user-specified parametric models for the estimands. For binary and categorical treatments, the model used for all estimands is of the form `E[Y|A=a,W] - E[Y|A=0,W] = formula(W)` where `formula` is specified by the user. For continuous treatments (only supported by `contglm`), the model used is `E[Y|A=a,W] - E[Y|A=0,W] = 1(a > 0) * formula_binary(W) + a * formula_continuous(W)`

causalglm also supports the following marginal structural model estimands:
 
1. Marginal structural models for the CATE: `E[CATE(W)|V] := E[E[Y|A=a,W] - E[Y|A=0,W]|V]` (categorical treatments)
2. Marginal structural models for the RR: `E[E[Y|A=a,W]|V]/E[E[Y|A=0,W]|V]` (categorical treatments)
3. Marginal structural models for the TSM : `E[E[Y|A=a,W]|V]` (categorical treatments)
4. Marginal structural models for the CATT : `E[CATE(W)|V, A=a] := E[E[Y|A=a,W] - E[Y|A=0,W]|V, A=a]` (categorical treatments)
 

causalglm consists of five main functions: 
 
1. `spglm` for semiparametric estimation of correctly specified parametric models for the `CATE`, `RR` and `OR`
2. `npglm` for robust nonparametric estimation of user-specified approximation models for the `CATE`, `CATT`, `TSM`, `RR` or `OR`
3. `msmglm` for robust nonparametric estimation of user-specified marginal structural models for the `CATE`, `CATT`, `TSM` or `RR`
4. `causalglmnet` for semiparametric estimation with high dimensional confounders `W` (a custom wrapper function for spglm focused on big data where standard ML may struggle)
5. `contglm` for robust nonparametric estimation of user-specified approximation models for the `CATE` as a function of a continuous or ordered numeric treatment.

The outputs of the methods include:

1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients
4. Individual-level treatment-effect predictions and 95\% confidence (prediction) intervals can be extracted with the `predict` function and argument `data`.
5. Plotting with `plot_msm` for objects returned by `msmglm`.

A rule of thumb for choosing between these methods is as follows:

1. Use `spglm` if you believe your parametric model for the treatment effect estimand is correct (this method is closest to glm)
2. Use `npglm` if you believe your parametric model for the treatment effect estimand is a good approximation but may not be correct (e.g. is missing some variables)
3. Use `msmglm` if you want to know how the treatment effect is causally affected by one or a number of variables `V` (fully adjusting for the remaining variables `W`) (or to learn univariate confounder-adjusted variable importance measures!)
4. Use `causalglmnet` if the variables `W` for which to adjust are (very) high dimensional.
5. Use `contglm` if your treatment is continuous or ordered and you are interested in the treatment effect per unit dose.

A longer answer is:

`spglm` is a semiparametric method which means that it assumes the user-specified parametric model is correct for inference. This method should be used if you are very confident in your parametric model. `npglm` is a nonparametric method that views the user-specified parametric model as an approximation or working-model for the true nonparametric estimand. The estimands are the best causal approximation of the true conditional estimand (i.e. projections). Because of this model agnostic view, npglm provides interpretable estimates and correct inference under no conditions. The user-specified parametric model need not be correct or even a good approximation for inference! `npglm` should be used if you believe your parametric model is a good approximation but are not very confident that it is correct. Also, it never hurts to run both `spglm` and `npglm` for robustness! If the parametric model is close to correct then the two methods should give similar estimates. Finally, `msmglm` deals with marginal structural models for the conditional treatment effect estimands. This method is useful if you are only interested in modeling the causal treatment effect as a function of a subset of variables `V` adjusting for all the available confounders `W` that remain. This allows for parsimonious causal modeling, still maximally adjusting for confounding. This function can be used to understand the causal variable importance of individual variables (by having `V` be a single variable) and allows for nice plots (see `plot_msm`).

 
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
 

### Conditional average treatment effect estimation


``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A * (1 + W + 2*W^2) + sin(5*W), sd = 0.3)
data <- data.frame(W,A,Y)

# Semiparametric
formula <- ~ poly(W, degree = 2, raw = FALSE)
output <-
  spglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "CATE",
    learning_method = "HAL",
    verbose = FALSE
  )

summary(output) 
head(predict(output, data = data))

# Nonparametric
formula <- ~ poly(W, degree = 2, raw = FALSE)
output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "CATE",
    learning_method = "HAL",
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
    learning_method = "HAL",
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
  W = "W", A = "A", Y = "Y"
)

summary(out)

# The CATE predictions are now a function of `A`
head(predict(out))

```


### Conditional odds ratio estimation
Note a log-linear model is used for the conditional odds ratio.
As a consequence, the parametric model specified by `formula` is actually a model for the log conditional odds ratio.

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
    estimand = "OR" 
  )
summary(output)

# Semiparametric inference
output <-
  spglm(
    ~1+W,
    data,
    W = c("W"), A = "A", Y = "Y",
    estimand = "OR" 
  )
summary(output)
```

 


### Conditional relative risk/treatment-effect estimation
 Note a log-linear model is used for the conditional relative risk.
As a consequence, the parametric model specified by `formula` is actually a model for the log conditional relative risk.

 

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

# Learn a marginal structural model with cool plots
output <-
  msmglm(
    formula,
    data,
    V = "W",
    W = "W", A = "A", Y = "Y",
    estimand = "RR",
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

 

## Learn marginal structural models for conditional treatment effects with `msmglm`

The `msmglm` function is a wrapper for `npglm` that focuses purely on marginal structural model estimation and it has some convenient plotting features. `V` can be multivariate but plotting is only supports for the univariate case. `estimand = "CATE"` corresponds with the estimand `E[CATE(W)|V]`, `estimand = "CATT"` corresponds with the estimand `E[CATE(W)|V, A=1]`, `estimand = "TSM"` corresponds with the estimand `E[E[Y|A=a,W]|V]`, and  `estimand = "RR"` corresponds with the estimand `E[E[Y|A=1,W]|V] / E[E[Y|A=0,W]|V]`.  The intercept model reduces to nonparametric efficient estimators for the marginal ATE, ATT, TSM and RR respectively. Note that formula should only depend on `V`. This method is useful to study confounder-adjusted associations between a continuous variable (i.e. `V`) and the treatment effect or outcome `Y`. 


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
 
 

