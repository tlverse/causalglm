# causalglm : interpretable and robust causal inference for heterogeneous treatment effects

## Installing
NOTE: This package is actively in development and is subject to continuous change. If you are unable to install it due to errors, try again in a day or two. Also, feel free to contact me personally or through the Issues tab.

 
 


To install this package, install the devtools CRAN package and run:

``` r
if(!require(devtools)) {
  install.packages(devtools)
}
devtools::install_github("tlverse/causalglm")
```

This package fully utilizes the powerful tlverse/tmle3 generalized targeted learning framework as well as the machine-learning frameworks tlverse/sl3 and tlverse/hal9001. This package also requires the following github packages which should be installed automatically.
``` r
devtools::install_github("tlverse/hal9001@devel")
devtools::install_github("tlverse/tmle3@general_submodels_devel")
devtools::install_github("tlverse/sl3@Larsvanderlaan-formula_fix")
```


## What is causalglm?

It is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for answers to causal questions, assuming parametric models can be dangerous. With even a seemingly small amount of confounding and misspecificaton, they can give biased answers. One way of mitigating this challenge is to only parametrically model the feature of the data-generating distribution that you care about. That is, use a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.

causalglm is an R package for robust generalized linear models and interpretable causal inference for heterogeneous (or conditional) treatment effects. Specifically, causalglm very significantly relaxes the assumptions needed for useful causal estimates and correct inference by employing semi and nonparametric models and adaptive machine-learning through targeted maximum likelihood estimation (TMLE) (van der Laan, Rose, 2011), while still robustly utilizing user-specified parametric forms for the conditional estimands and thus allowing for interpretable inference. Because of this, `causalglm` methods can (and often do in real-world settings) significantly reduce bias relative to conventional fully parametric methods like `glm`, and the cost in variance/confidence-interval-width is negligible. As another consequence, causalglm can be used in high dimensional settings and provides valid inference even when adaptive variable selection is used with, for example, the LASSO. See the writeup causalglm.pdf for a more theoretical overview of the methods implemented in this package. Also see the vignette for an overview of all the functionalities of causalglm. 

Throughout the development of causalglm, we greatly focused on making this package and all its methods as easy to use and understand as possible. Our goal is that this package can be used by a wide audience including amateurs, practioners, statisticians, causal-inference experts and nonexperts.

The statistical data-structure used throughout this package is `O = (W,A,Y)` where `W` represents a random vector of baseline (pretreatment) covariates/confounders, `A` is a usually binary treatment assignment with values in `c(0,1)`, and `Y` is some outcome variable. For marginal structural models, we also consider a subvector `V \subset W` that represents a subset of baseline variables that are of interest.

The estimands supported by causalglm are

1. Conditional average treatment effect (CATE) for arbitrary outcomes: `E[Y|A=1,W] - E[Y|A=0,W]`
2. Conditional odds ratio (OR) for binary outcomes: `{P(Y=1|A=1,W)/P(Y=0|A=1,W)}{P(Y=1|A=0,W)/P(Y=0|A=0,W)}`
3. Conditional relative risk (RR) for binary, count or nonnegative outcomes: `E[Y|A=1,W]/E[Y|A=0,W]`
4. Conditional treatment-specific mean (TSM) : `E[Y|A=a,W]`
5. Conditional average treatment effect among the treated (CATT) : the best approximation of `E[Y|A=1,W] - E[Y|A=0,W]` based on a user-specified formula/parametric model among the treated (i.e. observations with `A=1`)
 

causalglm also supports the following marginal structural model estimands:
 
1. Marginal structural models for the CATE: `E[CATE(W)|V] := E[E[Y|A=1,W] - E[Y|A=0,W]|V]`
2. Marginal structural models for the RR: `E[E[Y|A=1,W]|V]/E[E[Y|A=0,W]|V]`
3. Marginal structural models for the TSM : `E[E[Y|A=a,W]|V]`
4. Marginal structural models for the CATT : `E[CATE(W)|V, A=1] := E[E[Y|A=1,W] - E[Y|A=0,W]|V, A=1]`
 

causalglm consists of four main functions: 
 
1. `spglm` for semiparametric estimation of correctly specified parametric models for the `CATE`, `RR` and `OR`
2. `npglm` for robust nonparametric estimation for user-specified approximation models for the `CATE`, `CATT`, `TSM`, `RR` or `OR`
3. `msmglm` for robust nonparametric estimation for user-specified marginal structural models for the `CATE`, `CATT`, `TSM` or `RR`
4. `causalglmnet` for high dimensional confounders `W` (a custom wrapper function for spglm focused on big data where standard ML may struggle)
 
The outputs of the methods include:

1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients
4. Individual-level treatment-effect predictions and 95\% confidence (prediction) intervals can be extracted with the `predict` function and argument `data`.
5. Plotting with `plot_msm` for objects returned by `msmglm`.

A rule of thumb for choosing between these methods is as follows:

1. Use `spglm` if you believe your parametric model for the treatment effect estimand is correct (this method is closest to glm)
2. Use `npglm` if you believe your parametric model for the treatment effect estimand is a good approximation but may not be correct (or is missing some variables)
3. Use `msmglm` if you want to know how the treatment effect is causally affected by a one or a number of variables `V` (fully adjusting for the remaining variables `W`) (or to learn univariate confounder-adjusted variable importance measures!)
4. Use `causalglmnet` if the variables `W` for which to adjust are high dimensional.

A longer answer is:

`spglm` is a semiparametric method which means that it assumes the user-specified parametric model is correct for inference. This method should be used if you are very confident in your parametric model. `npglm` is a nonparametric method that views the user-specified parametric model as an approximation or working-model for the true nonparametric estimand. The estimands are the best causal approximation of the true conditional estimand (i.e. projections). Because of this model agnostic view, npglm provides interpretable estimates and correct inference under no conditions. The user-specified parametric model need not be correct or even a good approximation for inference! `npglm` should be used if you believe your parametric model is a good approximation but are not very confident that it is correct. Also, it never hurts to run both `spglm` and `npglm` for robustness! If the parametric model is close to correct then the two methods should give similar estimates. Finally, `msmglm` deals with marginal structural models for the conditional treatment effect estimands. This method is useful if you are only interested in modeling the causal treatment effect as a function of a subset of variables `V` adjusting for all the available confounders `W` that remain. This allows for parsimonious causal modeling, still maximally adjusting for confounding. This function can be used to understand the causal variable importance of individual variables (by having `V` be a single variable) and allows for nice plots (see `plot_msm`).

 
### User-friendly interface
 A minimalistic yet still very flexible front-end function for all routines is provided through the `spglm/npglm/causalglmnet/msmglm` functions. Check out the vignette to see how to use it! The necessary arguments are: 
1. A formula object for the `CATE`, `OR`, or `RR` (also `TSM`, `CATT` for `npglm`)
2. A data.frame containing the data
3. Variable names: `W`, `A`, `Y` are character vectors that store the variable names for the baseline variables, treatment variable and outcome variable.
4. Choice of estimand: `"CATE"`, `"OR"`, `"RR"` (also `"TSM"`, `"CATT"` for `npglm`)
 
 

### Conditional average treatment effect and partially-linear least-squares regression (spglm)
`spglm` with `estimand == "CATE"` performs estimation in the so-called "partially linear regression model" which *only* assumes

`CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] ~ a user-specified parametric model.`

This is equivalent to the semiparametric linear regression model `E[Y|A,W] = A CATE(W) + E[Y|A=0,W]` where `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W]` has a user-specified parametric form and `E[Y|A=0,W]` is an unspecified nuisance function that is learned nonparametrically using, for instance, machine-learning. In other words, only the treatment interaction terms in the linear model for `E[Y|A,W]` are modeled parametrically.  Here is some different ways you can model the CATE:

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A + W, sd = 0.3)
data <- data.frame(W,A,Y)


formula <- ~ 1 + W
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

# For spglm, you can pass in previous output objects to reuse the machine-learning fits.
# This only works if the new formula is a sub formula of the original formula (only for spglm)
formula <- ~ 1
output <-
  spglm(
    formula,
    output,
    estimand = "CATE",
    verbose = FALSE
  )

summary(output) 
head(predict(output, data = data))

 
```

### Conditional odds ratio and partially-linear logistic regression (spglm)
When Y is binary, the adjusted causal odds ratio between A and Y may be of interest. Use the function `spglm` with `estimand = "OR"`. The model used is the so-called "partially-linear logistic regression model" which *only* assumes

`logOR(W) := log[ {P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)} ] ~ user-specified parametric model`.

That is, the user specifies a parametric model for the log odds between A and Y and nothing else is assumed known. 

This is equivalent to assuming the semiparametric logistic regression model

`P(Y=1|A,W) = expit{A*logOR(W) + logit(P(Y=1|A=0,W))}`

where `P(Y=1|A=0,W)` is unspecified and learned using machine-learning. In other words, only the treatment interaction terms in the logistic regression model for `P(Y=1|A,W)` are modeled parametrically.

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n, size =  1, prob = plogis(A + A * W + W + sin(5 * W)))
data <- data.frame(W, A, Y)

formula <- ~ 1 + W
output <-
  spglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "OR" ,
    learning_method = "HAL",
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))
```

 


### Conditional relative risk/treatment-effect and partially-linear log-linear regression (spglm)
 
 When Y is binary, a count, or more generally nonnegative, the relative risk of Y with respect to A can be estimated. Use the function `spglm` with `estimand = "RR"`.

The model used is the so-called "partially-linear log-linear/poisson regression model" which *only* assumes

`log RR(W) := log{E[Y|A=1,W] / E[Y|A=0,W]} ~ user-specified parametric model`.

That is, we only assume the user specified parametric model (at the log scale) for the relative risk of Y with respect to A.

This is equivalent to assuming the semiparametric log-linear regression model
`E[Y|A,W] = E[Y|A=0,W] exp(log RR(W)) = E[Y|A=0,W] RR(W)`,
where `log RR(W)` is parametric and `E[Y|A=0,W]` is the background/placebo outcome model which is unspecified and learned using machine-learning. In other words, only the treatment interaction terms in the log-linear regression model for `E[Y|A,W]` are modeled parametrically.


``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp(A + A * W + sin(5 * W)))
data <- data.frame(W, A, Y)

formula <- ~ 1 + W
output <-
  spglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "RR" ,
    learning_method = "HAL",
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))
```


## Robust nonparametric inference for generalized linear models with npglm: CATE, CATT, TSM, RR, and OR

Rather than assuming a semiparametric model, we can instead make no assumptions on any functional forms (that is, assume a nonparametric model) and instead use a parametric or semiparametric model as an approximate "working model". This allows for interpretable coefficient-based estimates and inference that are correct under no assumptions on the functional form of the estimand. 

This nonparametric view is implemented in the function `npglm`. The estimates obtained are for the best approximation of the true estimand in the parametric "working model". That is, the estimands are the coefficients of the causal projection of the true estimand onto the parametric working model, where the projection will be defined next.  Even when you believe the working model is correct, this function may still be of interest for robustness. 

 

### Robust nonparametric inference for the CATE, CATT, and TSM (conditional treatment effects with npglm)

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A + W, sd = 0.3)
data <- data.frame(W,A,Y)

formula <- ~ 1 + W
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

# npglm (and msmglm) can reuse fits across formulas and estimands (no conditions!). Pass output to the data argument.
output <-
  npglm(
    formula,
   output,
    estimand = "CATT",
    verbose = FALSE
  )

summary(output) 
head(predict(output, data = data))


output <-
  npglm(
    formula,
    output,
    estimand = "TSM",
    verbose = FALSE
  )
# returns a list of causalglm objects for each level of `A`
output1 <- output[[1]]
output2 <- output[[2]]
summary(output1) 
head(predict(output1, data = data))
```
 

### Robust nonparametric inference for the OR (conditional odds ratio with npglm)

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n, size =  1, prob = plogis(A + A * W + W + sin(5 * W)))
data <- data.frame(W, A, Y)

formula <- ~ 1 + W
output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "OR" ,
    learning_method = "HAL",
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))
```


### Robust nonparametric inference for the RR (conditional relative-risk with npglm)
``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp(A + A * W + sin(5 * W)))
data <- data.frame(W, A, Y)

formula <- ~ 1 + W
output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "RR" ,
    learning_method = "HAL",
    verbose = FALSE
  )

summary(output)
head(predict(output, data = data))
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
 
 

