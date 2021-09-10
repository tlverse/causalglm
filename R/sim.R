


#' Simulate a dataset with known constant CATE
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.CATE <- function(n = 1500, p = 2, prop_active = 1, sigma = NULL, formula_estimand = ~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -1, 1)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A, as.data.frame(W))
  XY <- model.matrix(formula_Y0W, as.data.frame(W))
  V <- model.matrix(formula_estimand, as.data.frame(W))
  pA <- ncol(XA)
  pY <- ncol(XY)

  if (is.null(beta_A)) {
    activeA <- rbinom(pA, size = 1, prob = prop_active)
    beta_A <- runif(pA, min = -1, max = 1)
    beta_A <- 1.5 * beta_A * activeA / sum(abs(beta_A * activeA))
    names(beta_A) <- colnames(XA)
  }
  if (is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob = prop_active)
    beta_Y <- runif(pY, min = -1, max = 1)
    beta_Y <- 2 * beta_Y * activeY / sum(abs(beta_Y * activeY))

    names(beta_Y) <- colnames(XY)
  }

  g1 <- plogis(XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)
  Q0 <- 1 + XY %*% beta_Y

  if (!is.null(beta)) {
    beta_CATE <- beta
  } else {
    beta_CATE <- runif(ncol(V), min = -1, max = 1)
    beta_CATE <- beta_CATE / sum(abs(beta_CATE))
  }

  CATE <- V %*% beta_CATE

  Q <- CATE * A + Q0
  if (is.null(sigma)) {
    sigma <- sd(Q) / 4
  }
  Y <- rnorm(n, mean = Q, sd = sigma)

  data <- data.frame(W, A = A, Y = Y, pA1 = g1, pY = Q, sd = sigma, CATE = CATE)
  return(list(descr = "Data simulated from parametric linear model with constant conditional average treatment effect", beta_CATE = beta_CATE, data = data, W = colnames(W), A = "A", Y = "Y", beta_A = beta_A, beta_Y0 = beta_Y, link = "logistic"))
}
#' Simulate a dataset with known constant RR
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.RR <- function(n = 1500, p = 2, prop_active = 1, formula_estimand = ~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -0.5, 0.5)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A, as.data.frame(W))
  XY <- model.matrix(formula_Y0W, as.data.frame(W))
  V <- model.matrix(formula_estimand, as.data.frame(W))
  pA <- ncol(XA)
  pY <- ncol(XY)
  if (is.null(beta_A)) {
    activeA <- rbinom(pA, size = 1, prob = prop_active)
    beta_A <- runif(pA, min = -1, max = 1)
    beta_A <- 3 * beta_A * activeA / sum(abs(beta_A * activeA))
    names(beta_A) <- colnames(XA)
  }
  if (is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob = prop_active)
    beta_Y <- runif(pY, min = 0, max = 1)
    beta_Y <- 0.2 + 3 * beta_Y * activeY / sum(abs(beta_Y * activeY))
    names(beta_Y) <- colnames(XY)
  }
  g1 <- plogis(XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)
  Q0 <- exp(XY %*% beta_Y)
  if (!is.null(beta)) {
    betalogRR <- beta
  } else {
    betalogRR <- runif(ncol(V), min = 0, max = 1)
    betalogRR <- 1.4 * betalogRR / sum(abs(betalogRR))
  }

  RR <- exp(V %*% betalogRR)

  Q <- (1 - A) * Q0 + A * Q0 * RR
  Y <- rpois(n, lambda = Q)

  data <- data.frame(W, A = A, Y = Y, pA1 = g1, pY = Q, pY0 = Q0, pY1 = Q0 * RR, RR = RR)
  return(list(descr = "Data simulated from parametric linear model with known relative risk", beta_logRR = betalogRR, data = data, W = colnames(W), A = "A", Y = "Y", beta_A = beta_A, beta_Y0 = beta_Y, link = "logistic"))
}
#' Simulate a dataset with known constant OR
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.OR <- function(n = 1500, p = 2, prop_active = 1, formula_estimand = ~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -1, 1)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A, as.data.frame(W))
  XY <- model.matrix(formula_Y0W, as.data.frame(W))
  V <- model.matrix(formula_estimand, as.data.frame(W))
  pA <- ncol(XA)
  pY <- ncol(XY)
  if (is.null(beta_A)) {
    activeA <- rbinom(pA, size = 1, prob = prop_active)
    beta_A <- runif(pA, min = -1, max = 1)
    beta_A <- beta_A * activeA / sum(abs(beta_A * activeA))
    names(beta_A) <- colnames(XA)
  }
  if (is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob = prop_active)
    beta_Y <- runif(pY, min = -1, max = 1)
    beta_Y <- 3 * beta_Y * activeY / sum(abs(beta_Y * activeY))
    names(beta_Y) <- colnames(XY)
  }
  g1 <- plogis(XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)



  if (!is.null(beta)) {
    betalogOR <- beta
  } else {
    if (ncol(V) > 1) {
      betalogOR <- runif(ncol(V), min = -1, max = 1)
      betalogOR <- 1.2 * betalogOR / sum(abs(betalogOR))
      betalogOR <- c(betalogOR)
    } else {
      betalogOR <- 1.5
    }
  }

  logOR <- V %*% betalogOR
  Q <- plogis((A * logOR + XY %*% beta_Y))
  Q1 <- plogis((1 * logOR + XY %*% beta_Y))
  Q0 <- plogis((XY %*% beta_Y))
  Y <- rbinom(n, size = 1, prob = Q)

  data <- data.frame(W, A = A, Y = Y, pA1 = g1, pY = Q, pY1 = Q1, pY0 = Q0, OR = exp(logOR))
  return(list(descr = "Data simulated from parametric linear model with known odds ratio", beta_logOR = betalogOR, data = data, W = colnames(W), A = "A", Y = "Y", beta_A = beta_A, beta_Y0 = beta_Y, link = "logistic"))
}
