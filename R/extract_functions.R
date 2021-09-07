#' @export
summary.causalGLM <- function(object) {
  print(object$coefs)
  return(invisible(object$coefs))
}

#' @export
coef.causalGLM<- function(object) {
  out <- (object$coefs)
  out
}


#' @export
predict.causalGLM <- function(object, Wnew ) {
  n <- object$n
  Wnew <- as.matrix(Wnew)
  formula <- object$formula
  linkinv <- object$linkinv
  
  estimates <- as.matrix(object$coefs)[,"coefs"]
  var_scaled <- object$var_mat 
  
  V_newer <- model.matrix(formula, data = as.data.frame(Wnew))
  est_grid <-V_newer%*%estimates
  
  
  
  se_grid <- apply(V_newer,1, function(m) {
    sqrt(sum(m * (var_scaled %*%m)))
  } )
  Zvalue <- abs(sqrt(n) * est_grid/se_grid)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)
  
  estimates_transformed <- linkinv(est_grid)
  ci <- cbind(est_grid - 1.96*se_grid/sqrt(n), est_grid + 1.96*se_grid/sqrt(n))
  ci_transformed <- linkinv(ci)
  preds_new <- cbind(V_newer,  est_grid,   estimates_transformed, se_grid/sqrt(n),    ci_transformed,  Zvalue,  pvalue)
  colnames(preds_new) <- c(colnames(V_newer), "linear predictor", "estimate",  "se/sqrt(n)",  "CI_left",  "CI_right", "Z-score", "p-value" )
  return(preds_new)
  
}













