# causalGLM (in development but can be used with caution)
## Semiparametric and nonparametric generalized linear models for causal inference using Targeted Maximum Likelihood Estimation in low and high dimensions

It is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for answers to causal questions, assuming parametric models can be dangerous. With even a seemingly small amount of confounding and misspecificaton, they can give biased answers. One way of mitigating this challenge is to instead assume a parametric model for only the feature of the data-generating distribution that you care about. That is, assume a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question? It is not worth the risk of being wrong.

In this package, we utilize targeted machine-learning to generalize the parametric generalized linear models commonly used for treatment effect estimation (e.g. the R package glm) to the world of semi and nonparametric models. There is little-to-no loss in precision/p-values/confidence-interval-widths with these semiparametric methods relative to parametric generalized linear models, but the bias reduction from these methods can be substantial! Simulations suggest that these methods can work well with small sample sizes. We employ ensemble machine-learning (Super-Learning) that adapts the aggressiveness of the ML algorithms with sample size, thereby allowing for robust and correct inference in a diverse range of settings. All methods utilize targeted maximum likelihood estimation (TMLE) (van der Laan, Rose, 2012).

This package supports (semiparametric and nonparametric versions of) the following point-treatment estimands:

1. Conditional average treatment effect (CATE). (Causal semiparametric linear regression)
2. Conditional odds ratio (OR) between two binary variables. (Causal semiparametric logistic regression)
3. Conditional relative risk (RR) for nonnegative outcomes and a binary treatment. (Causal semiparametric log-linear relative-risk regression)
4. Conditional treatment-specific mean (TSM) for categorical treatments. (Only supported nonparametrically with causalRobustGLM)
5. Conditional average treatment effect among the treated (Only supported nonparametrically with causalRobustGLM)
6. Using causalRobustGLM with lower dimensional formula arguments, you can also learn marginal structural models for the CATE, CATT and RR.

This package also supports the following survival estimands:
1. Conditional hazard ratio between two treatments with `causalRobustCOXph`.


Each estimand can be modeled with a user-specified parametric model that is either assumed correct (`causalGLM` and `causalGLMnet`) or as an approximation, i.e. working model, of the nonparametric true estimand (`causalRobustGLM`). The former approach provides interpretable estimates and correct inference only when the parametric model is correct, and the latter approach provides interpretable estimates and nonparametrically correct inference even when the parametric model is incorrect.

Noticable features supported:
1. All methods utilize the powerful tlverse/tmle3 generalized targeted learning framework.
2. Interpretable semiparametric estimates and efficient inference even with adaptive estimation and variable selection.
3. High dimensional covariates and variable selection for confounders (with the wrapper function `causalGLMnet`).
4. General machine-learning tools with the tlverse/sl3 generalized machine-learning ecosystem.
5. Built-in machine-learning routines for diverse settings and immediate use.

### User-friendly interface
The functions are designed to be easy to use (any feedback will be greatly appreciated). A minimalistic yet still very flexible front-end function for all routines is provided through the `causalGLM/causalRobustGLM/causalGLMnet` functions. Check out the vignette to see how to use it! The necessary arguments are: 
1. A formula object for the `CATE`, `OR`, or `RR` (or `TSM`, `CATT` for `causalRobustGLM`)
2. A data.frame containing the data
3. Variable names: `W`, `A`, `Y` are character vectors that store the variable names for the baseline variables, treatment variable and outcome variable.
4. Choice of estimand: `"CATE"`, `"OR"`, `"RR"` (or `"TSM"`, `"CATT"` for `causalRobustGLM`)

That's it! Feel free to customize the machine-learning routines available using the "learning_method" argument. Built in options are: SuperLearner, HAL, glm, glmnet, gam, earth (MARS), CV-autotuned-xgboost. Cross-fitting is performed automatically. If you want to make your own learner, use the sl3_Learner argument and the tlverse/sl3 package.

Outputs include:
1. Coefficient estimates (using the S3 summary function)
2. Z-scores and p-values for coefficients (Still to come)
3. 95% confidence intervals for coefficients

## Semiparametric inference for generalized linear models with causalGLM: CATE, OR, and RR

The main function `causalGLM` implements semiparametric estimators for the CATE, OR and RR, which are each identified by some partially-linear generalized-linear model. We will utilize the statistical data-structure `O=(W,A,Y)` where `W` represents a vector of baseline variables, `A` is a binary treatment variable, and `Y` is some outcome.

### Conditional average treatment effect and partially-linear least-squares regression
`causalGLM` with `estimand == "CATE"` performs estimation in the so-called "partially linear regression model" defined as
`E[Y|A,W] = A CATE(W) + E[Y|A=0,W]` where `CATE(W) = E[Y|A=1,W] - E[Y|A=0,W]` has a user-specified parametric form and E[Y|A=0,W] is a nuisance function that is learned nonparametrically using machine-learning. Using the `formula` argument of `causalGLM`, one can learn the following CATE models:

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

### Conditional relative risk/treatment-effect and partially-linear log-linear regression
When Y is binary, a count, or more generally nonnegative, the relative risk of Y with respect to A can be estimated. Use the function `causalGLM` with `estimand = "RR"`.

The model used is the so-called "partially-linear log-linear/poisson regression model" which *only* assumes

`log RR(W) := log{E[Y|A=1,W] / E[Y|A=0,W]} ~ user-specified parametric model`.

That is, we only assume the user specified parametric model (at the log scale) for the relative risk of Y with respect to A.

This is equivalent to assuming the log-linear regression model
`E[Y|A,W] = E[Y|A=0,W] exp(log RR(W)) = E[Y|A=0,W] RR(W)`,
where log RR(W) is parametric and E[Y|A=0,W] is the background/placebo outcome model which is unspecified and learned using machine-learning.



## Robust nonparametric inference for generalized linear models with causalRobustGLM: CATE, CATT, TSM, RR, and OR
Rather than assuming a semiparametric model, we can instead make no assumptions (that is, assume a nonparametric model) and instead use a parametric or semiparametric model as an approximate "working model". This allows for interpretable coefficient-based estimates and inference that are correct under no assumptions on the functional form of the estimand. 

This nonparametric view is implemented in the function `causalRobustGLM`. The estimates obtained are for the best approximation of the true estimand in the parametric "working model". That is, the estimand are the coefficients of the projection of the true estimand onto the parametric working model, where the projection will be defined next.  Even when you believe the working model is correct, this function may still be of interest for robustness. For the most part, you can interpret the estimates in the same way you interpret the estimates given by `causalGLM`. 

We critically note that the semiparametric estimates given by `causalGLM` are (usually) not asymptotically equivalent to those given by`causalRobustGLM` when the parametric model is incorrect. The latter method can truly be viewed as an estimator for the best causal approximation, while the former is not necessarily so. Therefore, `causalRobustGLM` does not only give nonparametrically correct inference but also provides estimates for a nonparametric estimand that is more interpretable than the coefficients of the misspecified semiparametric model given by `causalGLM`. Notably, the intercept model for causalRobustGLM often corresponds with estimation of a nonparametric marginal causal parameter (like the ATE, ATT, marginal TSM, or marginal relative risk). This feature generalizes to marginal structural models for a number of the estimands. This is not true for the semiparametric methods implemented in `causalGLM`. There is a usually slight increase in confidence interval width for the nonparametric methods relative to the semiparametric methods.


### Robust nonparametric inference for the CATE
This method is useful for assessing heterogenity in the additive treatment effect across all individuals.

Let V := V(W) be the random vector obtained by applying the user-specified formula mapping to W (i.e. V = model.matrix(formula, W)). 

Consider the oracle least-squares risk function:

`R(beta) = E[(CATE(W) - beta^T * V )^2]`.

Our estimand of interest beta' is defined as the minimizer of the above risk function. In other words, we are interested in the coefficients of the least-squares projection onto the linear parametric working-model.

In particular, if V = 1 (i.e. `formula = ~1`) then the solution is equal to `beta := E[CATE(W)] = E[E[Y|A=1,W] - E[Y|A=0,W]] = AT`E, which is exactly the average treatment effect. Thus, this nonparametric working-model-based estimator does capture nonparametric causal additive treatment effects.

Notably, if `formula = ~1` is passed to `causalRobustGLM` then the coefficient is an efficient nonparametric estimator of the ATE, which may be of independent interest.

By specifying a formula of a lower dimensional feature `Z` of `W`, marginal structural models for the CATE can also be learned with this function. Specifically, if V(Z) is the design matrix obtained from the formula and one assumes
`E[CATE(W)|Z] = beta^T V(Z)`
then robustCausalGLM wil actually return estimates of the above beta (if the model is incorrect it can still be viewed as a working model approximation).

### Robust nonparametric inference for the CATT
This method is useful for assessing heterogenity in the additive treatment effect among the treated.

Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the oracle least-squares risk function:

`R(beta) = E(E[Y|A,W] - A * beta^T * V - E[Y|A=0,W])^2`,

where beta is a candidate coefficient vector for the best approximation of the true CATE. The minimizer beta' of the risk function `R(beta)` is our desired estimand.
This is equivalent to minimizing the risk function:

`R'(beta) = E[A * (CATE(W) - beta^T * V)]^2`,

which is the least-squares projection of `CATE(W) := E[Y|A=1,W] - E[Y|A=0,W]` onto the parametric working model `beta^T * V` using only observations with `A = 1` (among the treated). In particular, if `V = 1` (i.e. `formula = ~1`) then the solution is equal to `beta := E[CATE(W)|A=1] = E[E[Y|A=1,W] - E[Y|A=0,W]|A=1] = ATT`, which is exactly the average treatment effect among the treated. For this reason, we call this estimand the conditional average treatment-effect among the treated (CATT), since it is the best working-model approximation/predictor of the true CATE among the treated. This general working-model-based estimand can still be interpreted as a measure for the CATE, and be be preferred over the CATE method when there are positivity issues.

Notably, if `formula = ~1` is passed to `causalRobustGLM` then the coefficient is an efficient nonparametric estimator of the ATT, which may be of independent interest.

By specifying a formula of a lower dimensional feature `Z` of `W`, marginal structural models for the CATT can also be learned with this function. Specifically, if V(Z) is the design matrix obtained from the formula and one assumes
`E[CATE(W)|Z, A=1] = beta^T V(Z)`
then robustCausalGLM wil actually return estimates of the above beta (if the model is incorrect it can still be viewed as a working model approximation). 

### Robust nonparametric inference for the conditional TSM
This method is useful for assessing heterogeniety in the average outcome across all individuals for a given treatment intervention.

Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the oracle least-squares risk function:

`R(beta) = E(E[Y|A=a,W] -  beta^T * V )^2`,

Our estimand of interest is the risk minimizer, which is the least-squares projection of `TSM(W) := E[Y|A=a,W]` onto the parametric working model `beta^T * V`.  


Notably, if `formula = ~1` is passed to `causalRobustGLM` then the coefficient is an efficient nonparametric estimator of the marginal treatment specific mean `E_WE[Y|A=a,W]`, which may be of independent interest.  

### Robust nonparametric inference for the OR
This method is useful for assessing heterogeniety in the odds ratio.

Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the logistic working submodel:

`P_approx(Y=1|A,W) = expit{A * beta^T * V + logit(P(Y=1|A=0,W))}`.

Our estimand of interest `beta'` corresponds with the coefficient vector of the log-likelihood projection of the true distribution `P(Y=1|A,W)` onto the working submodel `P_approx(Y=1|A,W)`.


### Robust nonparametric inference for the RR
This method is useful for assessing heterogenity in the relative treatment effect.

Let V be the random vector obtained by applying the user-specified formula mapping to W. 

Consider the poisson log likelihood type risk function:

`R(beta) := E{E[Y|A=0,W] exp(beta^T V) - E[Y|A=1,W] beta^T V}`.

Our estimand of interest beta' corresponds with risk minimizer of the above risk function, which can be viewed as a log-linear projection of the relative risk onto the working model.

Notably, if formula = ~1 is passed to `causalRobustGLM` then the coefficient is an efficient nonparametric estimator of the log of the marginal relative risk, which may be of independent interest. That is, the estimand is exactly the log of `E_W E[Y|A=1,W] / E_W E[Y|A=0,W]`. 

More generally, this method can be used to learn marginal structural model parameters. Specifically, if one assumes the marginal structural model `log(E[E[Y|A=1,W]|Z]/E[E[Y|A=1,W]|Z]) = beta^T V(Z)` where `Z` is a subset of `W` and `V` is obtained from a formula that only depends on `Z`, then the coefficients can be learned by applying `causalRobustGLM` with the formula that gave `V(Z)`. This is true because, by conditioning, the risk function can be rewritten as
`R(beta)  = E{E[E[Y|A=0,W]|Z] exp(beta^T V(Z)) - E[E[Y|A=1,W]|Z] beta^T V(Z)}`.

## Interpretable nonparametric inference for the conditional hazard ratio between two treatments for censored time-to-event data

The function `causalRobustCOXph` allows you to estimate the parameters of a user-specified (time-dependent) parametric working-model for the conditional hazard ratio between two treatments. Specifically, this estimator is totally nonparametric and utilizes the approximate working model:

`log[P(T=t| T \geq t, A=1, W)/P(T=t| T \geq t, A=0, W)] = beta^T f(W,t)`

where f(W,t) is a user-specified vector-valued parametric function that depends on `W` and the time `t`. The coefficient vector `beta` is estimated by projecting a nonparametric estimator of the conditional hazard `P(T=t| T \geq t, A, W)` onto the exponential working model

`P_working(T=t| T \geq t, A, W) := exp(beta^T f(W,t)) P(T=t| T \geq t, A=0, W)`. 

This working model is the proportional hazards model except only for the treatment interaction term. There are no restrictions made on the placebo arm conditional hazard `P(T=t| T \geq t, A=0, W)`. The projection is obtained by minimizing the log-linear risk function:

`E[sum_{t=1,t_0} P(T=t| T \geq t, A=0, W) * exp{beta^T f(W,t)} - P(T=t| T \geq t, A=1, W) * beta^T f(W,t)]`,

which is essentially the same risk function used in the previous section for robust estimation of the relative risk.

By specifying lower dimensional formulas, marginal structural models for the hazard ratio can also be learned.

## Semiparametric inference for high dimensional generalized linear models with causalGLMnet (the LASSO): CATE, OR, and RR
For high dimensional W, you can use the wrapper function `causalGLMnet` which runs `causalGLM` using a custom glmnet-LASSO learner for estimation. This allows for robust and fast estimation in high dimensional settings where conventional machine-learning algorithms may struggle. Cross-fitting can be performed to reduce bias. This method can be viewed as an adaptive version of "glm" in that confounders/variables to adjust for are adaptively selected using the LASSO, while still allowing for asymptotically correct post-selection inference. 


## Need a new or specialized method? Questions? Suggestions?

Any confusion? Questions? Don't know which method to use? None of the methods handle your problem? Need a custom/specialized method?

Just send me a message. I would enjoy connecting with you and be happy to develop and implement a new method to add to this package.

## Future features to add

1. Categorical treatments (already supported by TSM)
2. Outcome missingness support
3. General link functions (mainly the poisson-link for the conditional TSM)
4. P-values and prediction confidence intervals.
5. Instrumental variable parameters
6. Survival parameters

If any of the above features is urgently needed (or just needed), let me know.

## References:
To be completed. 


These semiparametric models have a rich history and their theory goes back a long time. These references are very incomplete and and further references will be added in the future.

Most methods are based on theory and pseudo-code provided in the working paper van der Laan (2009), some of which is also published in journals: https://core.ac.uk/download/pdf/61320177.pdf (Page 600, 621, ish)


The relative risk method is treated in Targeted Maximum Likelihood Estimation of Conditional Relative Risk in a Semi-parametric Regression Model, Tuglus et al. (2011): https://biostats.bepress.com/ucbbiostat/paper283/.
The CATE method is treated in Statistical Inference for Variable Importance, van der Laan (2006): https://biostats.bepress.com/ucbbiostat/paper188/
For machine-learnng, the package tlverse/hal9001 and tlverse/sl3 are used: https://github.com/tlverse/hal9001

See also:

Estimation of a non-parametric variable importance measure of a continuous exposure, Chambaz et al. (2012): https://projecteuclid.org/journals/electronic-journal-Nonparametricof-statistics/volume-6/issue-none/Estimation-of-a-non-parametric-variable-importance-measure-of-a/10.1214/12-EJS703.full

Causal effects based on marginal structural models, Neugebauer, van der Laan (2007): 
https://www.researchgate.net/publication/222318646_Nonparametric_causal_effects_based_on_marginal_structural_models  

Related R packages: 

https://github.com/ck37/varimpact/tree/master/R
https://academic.oup.com/bioinformatics/article/31/18/3054/241218
https://cran.case.edu/web/packages/tmle.npvi/tmle.npvi.pdf

For fully nonparametric ATE-type methods, see the tmle package: https://cran.r-project.org/web/packages/tmle/index.html
Or tlverse/tmle3: https://tlverse.org


 
