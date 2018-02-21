
my_pcls_gam_R <- function(formula, data, lambda,  family, min_increase = 0, min_mu_k = 1e-16, type = "FlexGAM2") {
    
    if(type == "FlexGAM2n" && family$family == "gaussian") stop("FlexGAM2n not possible mit gaussian data. Use FlexGAM1n instead.")
    n <- nrow(data)
    formula_temp <- as.formula(formula[[3]])
    sm  <- smoothCon(formula_temp, data, knots=NULL)[[1]]
    
    Ain1 <- (-diag(ncol(sm$X)))[-ncol(sm$X),]
    Ain2 <- cbind(0,diag(ncol(sm$X)-1))
    AinA <- Ain1+Ain2
    if(type == "FlexGAM2"){
        if(family$family == "binomial") {
            Ain <- rbind(AinA,diag(ncol(sm$X)),-diag(ncol(sm$X)))
            bin <- c(rep(min_increase,times=nrow(AinA)),rep(0,times=ncol(sm$X)),rep(-1,times=ncol(sm$X)))
        }
        if(family$family == "poisson") {
            Ain <- rbind(AinA,diag(ncol(sm$X)))
            bin <- c(rep(min_increase,times=nrow(AinA)),rep(0,times=ncol(sm$X)))
        }
        if(family$family == "gaussian") {
            Ain <- rbind(AinA)
            bin <- c(rep(min_increase,times=nrow(AinA)))
        }
        if(family$family == "Gamma") {
            Ain <- rbind(AinA,diag(ncol(sm$X)))
            bin <- c(rep(min_increase,times=nrow(AinA)),rep(0,times=ncol(sm$X)))
        }
    }
    if(type == "FlexGAM2n"){
        if(family$family == "binomial") {
            Ain <- rbind(diag(ncol(sm$X)),-diag(ncol(sm$X)))
            bin <- c(rep(0,times=ncol(sm$X)),rep(-1,times=ncol(sm$X)))
        }
        if(family$family == "poisson") {
            Ain <- rbind(diag(ncol(sm$X)))
            bin <- c(rep(0,times=ncol(sm$X)))
        }
        if(family$family == "gaussian") {
            Ain <- rbind(NULL)
            bin <- c(rep(min_increase,times=nrow(NULL)))
        }
        if(family$family == "Gamma") {
            Ain <- rbind(diag(ncol(sm$X)))
            bin <- c(rep(0,times=ncol(sm$X)))
        }
    }
    
    M <- list()
    M$y <- data$y
    M$w <- data$y
    M$w[] <- 1
    M$X <- sm$X
    M$C <- matrix(0,0,0)
    M$S <- sm$S
    M$off <- 0
    M$sp <- lambda
    M$p <- seq(0.1,0.9,length=ncol(sm$X))
    M$Ain <- Ain
    M$bin <- bin
    
    check_condition <- M$Ain %*% M$p
    
    if(!all(check_condition[,1] >= M$bin)) {
        stop(paste0("Check condition \n \r  Vector of conditions: \n\r",paste(check_condition[,1] >= M$bin, collapse=" "))) 
        }
    
    beta_1 <- pcls(M)
    
    fitted <- M$X %*% beta_1
    
    list("M"=M, "coefficients"=beta_1, "fitted"=fitted[,1], "formula" = formula, "sm" = sm, "type" = type)
}









my_pcls_gam <- function(formula, data, lambda = 10, family, min_increase=0, min_mu_k = 1e-16, type = "FlexGAM2") {
    
    bgf <- my_pcls_gam_R(formula = formula, data = data, lambda = lambda, family=family, min_increase = min_increase, min_mu_k = min_mu_k, type = type)
    list("coefficients"=bgf$coef, "fitted"=bgf$fitted, "linear.predictor" = bgf$fitted, "sp"=lambda, 
         "score"=bgf$score, "bgf"=bgf, "formula" = formula, "family" = family, "type" = type, "eta_k" = data$eta_k)
}


predict_my_pcls_gam <- function(f_k, data, type = "FlexGAM2") {
    eta_k <- data$eta_k
    X_new <- splines::splineDesign(f_k$bgf$sm$knots, eta_k, ord = 4, 
                                   rep(0,times=length(eta_k)), outer.ok = TRUE,
                                    sparse = FALSE)
    fitted <- (X_new %*% f_k$coef)[,1]
    fitted[eta_k < min(f_k$eta_k)] <- f_k$fitted[which.min(f_k$eta_k)]
    fitted[eta_k > max(f_k$eta_k)] <- f_k$fitted[which.max(f_k$eta_k)]
    fitted
}
