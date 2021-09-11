# causalglm 

To install this package, install the devtools CRAN package and run:

``` r
devtools::install_github("tlverse/causalglm")
```

For a in-depth description of these methods and example code, see the document "causalglm_writeup.pdf" in the "writeup" folder. This readme is largely a condensed version of this writeup. For example code and a walk-through guide, also see the "vignette.Rmd" document in the "vignette" folder.  

## Semiparametric and nonparametric generalized linear models for conditional causal inference using Targeted Maximum Likelihood Estimation

  
It is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for answers to causal questions, assuming parametric models can be dangerous. With even a seemingly small amount of confounding and misspecificaton, they can give biased answers. One way of mitigating this challenge is to instead assume a parametric model for only the feature of the data-generating distribution that you care about. That is, assume a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.

In this package, we utilize targeted machine-learning to generalize the parametric generalized linear models commonly used for treatment effect estimation (e.g. the R package glm) to the world of semi and nonparametric models. There is little-to-no loss in precision/p-values/confidence-interval-widths with these semiparametric methods relative to parametric generalized linear models, but the bias reduction from these methods can be substantial! Simulations suggest that these methods can work well with small sample sizes. We employ ensemble machine-learning (Super-Learning) that adapts the aggressiveness of the ML algorithms with sample size, thereby allowing for robust and correct inference in a diverse range of settings. All methods utilize targeted maximum likelihood estimation (TMLE) (van der Laan, Rose, 2012).

This package supports (semiparametric and nonparametric versions of) the following point-treatment estimands:

1. Conditional average treatment effect (CATE). (Causal semiparametric linear regression)
2. Conditional odds ratio (OR) between two binary variables. (Causal semiparametric logistic regression)
3. Conditional relative risk (RR) for nonnegative outcomes and a binary treatment. (Causal semiparametric log-linear relative-risk regression)
4. Conditional treatment-specific mean (TSM) for categorical treatments. (Only supported nonparametrically with npglm)
5. Conditional average treatment effect among the treated (CATT) (Only supported nonparametrically with npglm)
6. Using npglm with lower dimensional formula arguments, you can also learn marginal structural models for the CATE, CATT, TSM and RR.

 
The semiparametric methods are run using the function `spglm` and the nonparametric methods are run using the function `npglm`. 
A semiparametric high dimensional LASSO version of `spglm` is implemented in `causalglmnet`. 

Each estimand can be modeled with a user-specified parametric model that is either assumed correct (`spglm` and `causalglmnet`) or as an approximation, i.e. working model, of the nonparametric true estimand (`npglm`). The former approach provides interpretable estimates and correct inference only when the parametric model is correct, and the latter approach provides interpretable estimates and nonparametrically correct inference even when the parametric model is incorrect.

 
### User-friendly interface
The functions are designed to be easy to use (any feedback will be greatly appreciated). A minimalistic yet still very flexible front-end function for all routines is provided through the `spglm/npglm/causalglmnet` functions. Check out the vignette to see how to use it! The necessary arguments are: 
1. A formula object for the `CATE`, `OR`, or `RR` (also `TSM`, `CATT` for `npglm`)
2. A data.frame containing the data
3. Variable names: `W`, `A`, `Y` are character vectors that store the variable names for the baseline variables, treatment variable and outcome variable.
4. Choice of estimand: `"CATE"`, `"OR"`, `"RR"` (also `"TSM"`, `"CATT"` for `npglm`)

That's it! Feel free to customize the machine-learning routines available using the "learning_method" argument. Built in options are: SuperLearner, HAL, glm, glmnet, gam, earth (MARS), CV-autotuned-xgboost. Cross-fitting is performed automatically. If you want to make your own learner, use the sl3_Learner argument and the tlverse/sl3 package.

Outputs include:
1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients

## Semiparametric inference for generalized linear models with spglm: CATE, OR, and RR  

The function `spglm` implements semiparametric estimators for the CATE, OR and RR, which are each identified by some partially-linear generalized-linear model. We will utilize the statistical data-structure `O=(W,A,Y)` where `W` represents a vector of baseline variables, `A` is a binary treatment variable, and `Y` is some outcome.

### Conditional average treatment effect and partially-linear least-squares regression (spglm)
`spglm` with `estimand == "CATE"` performs estimation in the so-called "partially linear regression model" defined as
`E[Y|A,W] = A CATE(W) + E[Y|A=0,W]` where `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W]` has a user-specified parametric form and `E[Y|A=0,W]` is a nuisance function that is learned nonparametrically using machine-learning. Here is some different ways you can model the CATE:

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A + W, sd = 0.3)
data <- data.frame(W,A,Y)

formula <- ~ 1
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
```

### Conditional odds ratio and partially-linear logistic regression (spglm)
When Y is binary, the adjusted causal odds ratio between A and Y may be of interest. Use the function `spglm` with `estimand = "OR"`. The model used is the so-called "partially-linear logistic regression model" which *only* assumes

`logOR(W) := log[ {P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)} ] ~ user-specified parametric model`.

That is, the user specifies a parametric model for the log odds between A and Y and nothing else is assumed known. 

This is equivalent to assuming the logistic regression model

`P(Y=1|A,W) = expit{A*logOR(W) + logit(P(Y=1|A=0,W))}`

where `P(Y=1|A=0,W)` is unspecified and learned using machine-learning.

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n, size =  1, prob = plogis(A + A * W + W + sin(5 * W)))
data <- data.frame(W, A, Y)
formula ~ 1 + W
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
```

 


### Conditional relative risk/treatment-effect and partially-linear log-linear regression (spglm)
 
 When Y is binary, a count, or more generally nonnegative, the relative risk of Y with respect to A can be estimated. Use the function `spglm` with `estimand = "RR"`.

The model used is the so-called "partially-linear log-linear/poisson regression model" which *only* assumes

`log RR(W) := log{E[Y|A=1,W] / E[Y|A=0,W]} ~ user-specified parametric model`.

That is, we only assume the user specified parametric model (at the log scale) for the relative risk of Y with respect to A.

This is equivalent to assuming the log-linear regression model
`E[Y|A,W] = E[Y|A=0,W] exp(log RR(W)) = E[Y|A=0,W] RR(W)`,
where `log RR(W)` is parametric and `E[Y|A=0,W]` is the background/placebo outcome model which is unspecified and learned using machine-learning.


``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp(A + A * W + sin(5 * W)))
data <- data.frame(W, A, Y)
formula ~ 1 + W
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
```


## Robust nonparametric inference for generalized linear models with npglm: CATE, CATT, TSM, RR, and OR
Rather than assuming a semiparametric model, we can instead make no assumptions (that is, assume a nonparametric model) and instead use a parametric or semiparametric model as an approximate "working model". This allows for interpretable coefficient-based estimates and inference that are correct under no assumptions on the functional form of the estimand. 

This nonparametric view is implemented in the function `npglm`. The estimates obtained are for the best approximation of the true estimand in the parametric "working model". That is, the estimands are the coefficients of the projection of the true estimand onto the parametric working model, where the projection will be defined next.  Even when you believe the working model is correct, this function may still be of interest for robustness. 

We critically note that the semiparametric estimates given by `spglm` are (usually) not asymptotically equivalent to those given by`npglm` when the parametric model is incorrect. The latter method can truly be viewed as an estimator for the best causal approximation, while the former is not necessarily so. Therefore, `npglm` does not only give nonparametrically correct inference but also provides estimates for a nonparametric estimand that is more interpretable than the coefficients of the misspecified semiparametric model given by `spglm`. Notably, the intercept model for npglm often corresponds with estimation of a nonparametric marginal causal parameter (like the ATE, ATT, marginal TSM, or marginal relative risk). This feature generalizes to marginal structural models for a number of the estimands. This is not true for the semiparametric methods implemented in `spglm`. There is a usually slight increase in confidence interval width for the nonparametric methods relative to the semiparametric methods.

We refer to the writeup "causalglm_writeup.pdf", which can be found in the folder "writeup" for the description of the following working-model based methods.

 

### Robust nonparametric inference for the CATE, CATT, and TSM (conditional treatment effects with npglm)

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rnorm(n, mean = A + W, sd = 0.3)
data <- data.frame(W,A,Y)

formula <- ~ 1
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

output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "CATT",
    learning_method = "HAL",
    verbose = FALSE
  )
summary(output) 


output <-
  npglm(
    formula,
    data,
    W = "W", A = "A", Y = "Y",
    estimand = "TSM",
    learning_method = "HAL",
    verbose = FALSE
  )
summary(output) 
```

### Robust nonparametric inference for the OR (conditional odds ratio with npglm)

``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n, size =  1, prob = plogis(A + A * W + W + sin(5 * W)))
data <- data.frame(W, A, Y)
formula ~ 1 + W
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
```


### Robust nonparametric inference for the RR (conditional relative-risk with npglm)
``` r
library(causalglm)
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp(A + A * W + sin(5 * W)))
data <- data.frame(W, A, Y)
formula ~ 1 + W
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

```

## Need a new or specialized method? Questions? Suggestions?

Any confusion? Questions? Don't know which method to use? None of the methods handle your problem? Need a custom/specialized method?
  
  Just send me a message. I would enjoy connecting with you and be happy to develop and implement a new method to add to this package. Please contact me at vanderlaanlars@yahoo.com.

 
## References (see writeup)
 
 

