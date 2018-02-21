coefficients.flexgam <- function(object, ... ) {
    coef <- c(object$beta_k,object$f_k$coefficients)
    names(coef) <- c(names(object$beta_k), paste("Psi.",1:length(object$f_k$coefficients),sep=""))
    coef
}

coef.flexgam <- function(object, ... ) {
    coef <- c(object$beta_k,object$f_k$coefficients)
    names(coef) <- c(names(object$beta_k), paste("Psi.",1:length(object$f_k$coefficients),sep=""))
    coef
}

residuals.flexgam <- function(object, ...) {
    res <- object$response - object$fitted
    attr(res,"names") <- NULL
    as.numeric(res)
}

resid.flexgam <- function(object, ...) {
    res <- object$response - object$fitted
    attr(res,"names") <- NULL
    as.numeric(res)
}

fitted.flexgam <- function(object, ...) {
    as.numeric(object$fitted)
}

fitted.values.flexgam <- function(object, ...) {
    as.numeric(object$fitted)
}
