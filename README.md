# causalGLM
## Semiparametric and nonparametric generalized linear models for causal inference using Targeted Maximum Likelihood Estimation 

It is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for answers to causal questions, assuming parametric models can be dangerous. With even a seemingly small amount of confounding and misspecificaton, they can give biased answers. One way of mitigating this challenge is to instead assume a parametric model for only the feature of the data-generating distribution that you care about. That is, assume a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.

In this package, we utilize targeted machine-learning to generalize the parametric generalized linear models commonly used for treatment effect estimation (e.g. the R package glm) to the world of semi and nonparametric models. There is little-to-no loss in precision/p-values/confidence-interval-widths with these semiparametric methods relative to parametric generalized linear models, but the bias reduction from these methods can be substantial! Simulations suggest that these methods can work well with small sample sizes. We employ ensemble machine-learning (Super-Learning) that adapts the aggressiveness of the ML algorithms with sample size, thereby allowing for robust and correct inference in a diverse range of settings. All methods utilize targeted maximum likelihood estimation (TMLE) (van der Laan, Rose, 2012).

This package supports (semiparametric and nonparametric versions of) the estimands:

1. Conditional average treatment effect (CATE) estimation. (Causal semiparametric linear regression)
2. Conditional odds ratio (OR) estimation between two binary variables]. (Causal semiparametric logistic regression)
3. Conditional relative risk (RR) regression for nonnegative outcomes and a binary treatment. (Causal semiparametric log-linear relative-risk regression)


Noticable features supported:
1. Interpretable semiparametric estimates and inference even with adaptive estimation and variable selection.
2. High dimensional covariates and variable selection for confounders (with the wrapper function causalGLMnet).
3. General machine-learning tools with the tlverse/sl3 ecosystem.
4. Built-in machine-learning routines for diverse settings and immediate use.

### User-friendly interface
The functions are designed to be easy to use (any feedback will be greatly appreciated). A minimalistic yet still very flexible front-end function for all routines is provided through the causalGLM/causalRobustGLM/causalGLMnet functions. Check out the vignette to see how to use it! The necessary arguments are: 
1. A formula object for the CATE, OR, or RR
2. A data.frame containing the data
3. Character vectors: W, A, Y that store the variable names for the baseline variables, treatment variable and outcome variable.
4. Choice of estimand: "CATE", "OR", "RR"

That's it! Feel free to customize the machine-learning routines available using the "learning_method" argument. Built in options are: SuperLearner, HAL, glm, glmnet, gam, earth (MARS), CV-autotuned-xgboost. Cross-fitting is performed automatically. If you want to make your own learner, use the sl3_Learner argument and the tlverse/sl3 package.

Outputs include:
1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients (Still to come)
3. 95% confidence intervals for coefficients

## Semiparametric inference for generalized linear models with causalGLM: CATE, OR, and RR

The main function `causalGLM` implements semiparametric estimators for the CATE, OR and RR, which are each identified by some partially-linear generalized-linear model. We will utilize the statistical data-structure `O=(W,A,Y)` where `W` represents a vector of baseline variables, `A` is a binary treatment variable, and `Y` is some outcome.

### Conditional average treatment effect and partially-linear least-squares regression
causalGLM with `estimand == "CATE"` performs estimation in the so-called "partially linear regression model" defined as
E[Y|A,W] = A CATE(W) + E[Y|A=0,W] where CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] has a user-specified parametric form and E[Y|A=0,W] is a nuisance function that is learned nonparametrically using machine-learning. Using the `formula` argument of causalGLM, one can learn the following CATE models:

1. Constant CATE: `formula = ~ 1`
This formula encodes the model `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] = a` for some constant coefficient `a`.

2. Effect modification and subgroup effects: `formula = ~ 1 + W_1`
This formula encodes the model `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] = a + b * W` for coefficients `a` and `b`.

3. More complex: `formula = ~ 1 + W_1 + W_1 * W_2`
This formula encodes the model `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] = a + b * W + c * W_1 * W_2` for coefficients `a`, `b` and `c`.


### Conditional odds ratio and partially-linear logistic regression
When Y is binary, the adjusted causal odds ratio between A and Y may be of interest. Use the function `causalGLM` with `estimand = "OR"`.


The model used is the so-called "partially-linear logistic regression model" which *only* assumes

`logOR(W) := log[ {P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)} ] ~ user-specified parametric model`.
That is, the user specifies a parametric model for the log odds between A and Y and nothing else is assumed known.

This is equivalent to assuming the logistic regression model

`P(Y=1|A,W) = expit{A*logOR(W) + logit(P(Y=1|A=0,W))}`

where P(Y=1|A=0,W) is unspecified and learned using machine-learning.

Just like with the CATE, you can specify arbitrary parametric forms of the conditional odds ratio (at the log-scale).

### Conditional relative risk/treatment-effect and partially-linear log-linear/link regression
When Y is binary, a count, or more generally nonnegative, the relative risk of Y with respect to A can be estimated. Use the function `causalGLM` with `estimand = "RR"`.

The model used is the so-called "partially-linear log-linear/poisson regression model" which *only* assumes

`log RR(W) := log{E[Y|A=1,W] / E[Y|A=0,W]} ~ user-specified parametric model`.

That is, we only assume the user specified parametric model (at the log scale) for the relative risk of Y with respect to A.

This is equivalent to assuming the log-linear regression model
`E[Y|A,W] = E[Y|A=0,W] exp(log RR(W)) = E[Y|A=0,W] RR(W)`,
where log RR(W) is parametric and E[Y|A=0,W] is the background/placebo outcome model which is unspecified and learned using machine-learning.



## Robust nonparametric inference for generalized linear models with causalRobustGLM: CATE, CATT, and OR
Rather than assuming a semiparametric model, we can instead make no assumptions (that is, assume a nonparametric model) and instead use a parametric or semiparametric model as an approximate "working model". This allows for interpretable coefficient-based estimates and inference that are correct under no assumptions on the functional form of the estimand. 

This nonparametric view is implemented in the function `causalRobustGLM`. The estimates obtained are for the best approximation of the true estimand in the parametric "working model". That is, the estimand are the coefficients of the projection of the true estimand onto the parametric working model, where the projection will be defined next.  Even when you believe the working model is correct, this function may still be of interest for robustness. For the most part, you can interpret the estimates in the same way you interpret the estimates given by `causalGLM`. We now will define the actual 

### Robust nonparametric inference for the CATE
Let V := V(W) be the random vector obtained by applying the user-specified formula mapping to W (i.e. V = model.matrix(formula, W)). 

Consider the oracle least-squares risk function:

`R(beta) = E[(CATE(W) - beta^T * V )^2]`.

Our estimand of interest beta' is defined as the minimizer of the above risk function. In other words, we are interested in the coefficients of the least-squares projection onto the linear parametric working-model.

In particular, if V = 1 (i.e. formula = ~1) then the solution is equal to `beta := E[CATE(W)] = E[E[Y|A=1,W] - E[Y|A=0,W]] = AT`E, which is exactly the average treatment effect. Thus, this nonparametric working-model-based estimator does capture nonparametric causal additive treatment effects.

Notably, if `formula = ~1` is passed to causalRobustGLM then the coefficient is an efficient nonparametric estimator of the ATE, which may be of independent interest.

### Robust nonparametric inference for the CATT
Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the oracle least-squares risk function:

`R(beta) = E(E[Y|A,W] - A * beta^T * V - E[Y|A=0,W])^2`,

where beta is a candidate coefficient vector for the best approximation of the true CATE. The minimizer beta' of the risk function R(beta) is our desired estimand.
This is equivalent to minimizing the risk function:

`R'(beta) = E[A * (CATE(W) - beta^T * V)]^2`,

which is the least-squares projection of `CATE(W) := E[Y|A=1,W] - E[Y|A=0,W]` onto the parametric working model beta^T * V using only observations with A = 1 (among the treated). In particular, if V = 1 (i.e. formula = ~1) then the solution is equal to `beta := E[CATE(W)|A=1] = E[E[Y|A=1,W] - E[Y|A=0,W]|A=1] = ATT`, which is exactly the average treatment effect among the treated. For this reason, we call this estimand the conditional average treatment-effect among the treated (CATT), since it is the best working-model approximation/predictor of the true CATE among the treated. This general working-model-based estimand can still be interpreted as a measure for the CATE, and be be preferred over the CATE method when there are positivity issues.

Notably, if formula = ~1 is passed to causalRobustGLM then the coefficient is an efficient nonparametric estimator of the ATT, which may be of independent interest.


### Robust nonparametric inference for the OR
Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the logistic working submodel:

`P_approx(Y=1|A,W) = expit{A * beta^T * V + logit(P(Y=1|A=0,W))}`.

Our estimand of interest beta' corresponds with the coefficient vector of the log-likelihood projection of the true distribution P(Y=1|A,W) onto the working submodel P_approx(Y=1|A,W).

## Semiparametric inference for high dimensional generalized linear models with causalGLMnet (the LASSO): CATE, OR, and RR
For high dimensional W, you can use the wrapper function "causalGLMnet" which runs "causalGLM" using a custom glmnet-LASSO learner for estimation. This allows for robust and fast estimation in high dimensional settings where conventional machine-learning algorithms may struggle. Cross-fitting can be performed to reduce bias. This method can be viewed as an adaptive version of "glm" in that confounders/variables to adjust for are adaptively selected using the LASSO, but still allow for asymptotically correct post-selection inference. 


## Need a new or specialized method?

Any confusion? Questions? Don't know which method to use? None of the methods handle your problem? Need a custom/specialized method?

Just send me a message. I would enjoy connecting with you and be happy to develop and implement a new method to add to this package.

