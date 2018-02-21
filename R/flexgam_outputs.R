
predict.flexgam <- function(object, newdata=NULL, type=c("response","linear.predictor","terms"), ...)  {
    
    type2 <- match.arg(type)
    
    norm_beta      <- FALSE               # object$control$norm_beta
    type           <- object$control$type
    gam_temp       <- object$gam_temp
    mean_eta_k     <- object$mean_eta_k
    sd_eta_k       <- object$sd_eta_k
    colMeans_f     <- object$colMeans_f
    sums_all_f     <- object$sums_all_f
    beta_k         <- object$beta_k
    f_k            <- object$f_k
    data_pred_norm <- object$data_pred_norm
    
    # rowSums(predict(object,type == "terms" )) != predict(object$linear.predictor) 
    # since linear.predictor is scaled!
    
    if(!is.null(newdata)) {
        if(!object$details$conv_out_first) {
            if(!norm_beta) {
                X_gam <- predict(gam_temp, newdata = newdata, type="lpmatrix")
                eta_k <- X_gam %*% beta_k
                eta_k <- ( eta_k - mean_eta_k) / sd_eta_k
                terms1 <- predict(gam_temp, newdata = newdata, type="terms")    
                
            } else {
                effects_per_f <- predict(gam_temp, newdata = newdata, type = "terms")
                
                
                # new_effects_per_f_temp <- (effects_per_f)/sums_all_f
                # new_effects_per_f  <- terms1 <- new_effects_per_f_temp #- colMeans_f
                # eta_k <-  rowSums(new_effects_per_f)
                # 
                
                
                
                sums_all_f <- sqrt(sum(predict(gam_temp, type = "terms", newdata=data_pred_norm)^2))
                
                new_effects_per_f_temp <- (effects_per_f)/sums_all_f
                colMeans_f <- colMeans(new_effects_per_f_temp)
                new_effects_per_f  <- new_effects_per_f_temp # - colMeans_f1
                
                eta_k <- rowSums(new_effects_per_f)
                
                mean_eta_k <- mean(eta_k)
                sd_eta_k   <- sd(eta_k)
                
                
                
            }
            
            if(type == "FlexGAM2" || type == "FlexGAM2n") {
                fitted <- (Psi_all(eta=eta_k, f_k=f_k, type=type, family = object$family))
            } 
            if(type == "FlexGAM1" || type == "FlexGAM1n") {
                fitted <-  Psi_all(eta=eta_k, f_k=f_k, type=type, family = object$family, g_inv = TRUE)
            }
        } else {
            fitted <- predict(object$gam_orig, newdata=newdata, type = "response")
            eta_k  <- predict(object$gam_orig, newdata=newdata, type = "link")
            terms1 <- predict(object$gam_orig, newdata=newdata, type = "terms")    
            
            if(!norm_beta) {
                eta_k <- ( eta_k - mean_eta_k) / sd_eta_k
            }
        }
    } else {
        fitted <- object$fitted
        eta_k  <- object$eta_k
        terms1 <- predict(gam_temp, type="terms")    
    }
    if(type2 == "response") {
        res <- as.numeric(fitted)
    } 
    if(type2 == "linear.predictor") {
        res <- as.numeric(eta_k)
    }
    if(type2 == "terms") {
        res <- terms1
    }
    res
}

deviance.flexgam <- function(object, newdata=NULL, use_penalty=FALSE, ...){
    if(is.null(use_penalty)) stop("TRUE/FALSE needed for use_penalty")
    if(is.na(use_penalty)) stop("TRUE/FALSE needed for use_penalty")
    
    if(is.null(newdata)) {
        fitted <- predict(object=object)
        y <- object$response
    } else {
        fitted <- predict(object=object, newdata=newdata)
        y <- newdata[,as.character(object$formula)[2]]
    }
    if(object$family$family == "binomial") {
        if((any(fitted >= 1-1e-16) || any(fitted <    1e-16)) ) warning("Predicted 0 or 1 occured and were slighty changed to avoid NAs") 
        if(any(fitted >= 1-1e-16)) fitted[fitted >= 1-1e-16] <- 1-1e-16 # Since predict.gam uses linear extrapolation it might occure 
        if(any(fitted <    1e-16)) fitted[fitted <    1e-16] <-   1e-16 # that the predicted values are outside the range of [0,1].
        dev <- -2*(sum(log(1-fitted[y==0])) + sum(log(fitted[y==1])))
    }
    if(object$family$family == "poisson") {
        if( any(fitted <    1e-16) ) warning("Predicted 0 occured and were slighty changed to avoid NAs") 
        if(any(fitted <    1e-16)) fitted[fitted <    1e-16] <-   1e-16 
        dev <- 2*(sum(y[y>0] * log(y[y>0] / fitted[y>0])) + sum( - (y - fitted)))
    }
    if(object$family$family == "gaussian") {
        dev <-  sum(  (y - fitted)^2)
    }
    if(object$family$family == "Gamma") {
        if( any(fitted <    1e-16)  ) warning("Predicted 0 occured and were slighty changed to avoid NAs") 
        if(any(fitted <    1e-16)) fitted[fitted <    1e-16] <-   1e-16 
        dev <- 2*(sum( (y - fitted) / fitted -log(y/fitted)))
    }
    
    if(use_penalty) {
        gam_penalty_vec <- 0
        if(length(object$gam_temp$smooth) > 0) {
            gam_penalty_vec <- NULL
            for(jj in 1: length(object$gam_temp$smooth)){
                gam_S_temp  <- object$gam_temp$smooth[[jj]]$S
                gam_sp_temp <- object$gam_temp$smooth[[jj]]$sp
                gam_S_temp1 <- gam_S_temp[[1]]
                gam_S_temp1[] <- 0
                for(jj_1 in 1:length(gam_S_temp)) {
                    gam_S_temp1 <- gam_S_temp1 + gam_S_temp[[jj_1]] * gam_sp_temp[jj_1]
                }
                if(any(gam_sp_temp == -1)) stop("Major error sp should be fixed!")
                gam_coef_temp <- object$beta_k[grepl(x=names(object$beta_k),pattern=object$gam_temp$smooth[[jj]]$label,fixed=TRUE)]
                #gam_coef_temp <- beta_k[1:19 + (19*(jj-1))]
                gam_penalty_vec <- c(gam_penalty_vec, (t(gam_coef_temp) %*% gam_S_temp1 %*% gam_coef_temp))
                }
            }
        f_k_penalty <- 0
        if(object$type == "FlexGAM1" || object$type == "FlexGAM1n"){
            f_k_S_temp  <- object$f_k$smooth[[1]]$S[[1]]
            f_k_sp_temp <- object$f_k$smooth[[1]]$sp
            if(f_k_sp_temp == -1) f_k_sp_temp <- object$f_k$sp[names(f_k_sp_temp)]
            f_k_coef_temp <- object$f_k$coefficients[grepl(x=names(object$f_k$coefficients),pattern=names(f_k_sp_temp),fixed=TRUE)]
            f_k_penalty <- f_k_sp_temp * (t(f_k_coef_temp) %*% f_k_S_temp %*% f_k_coef_temp)
        }
        if(object$type == "FlexGAM2" || object$type == "FlexGAM2n") {
            f_k_S_temp  <- object$f_k$bgf$sm$S[[1]]
            f_k_sp_temp <- object$f_k$sp
            f_k_coef_temp <- object$f_k$coefficients
            f_k_penalty <- f_k_sp_temp * (t(f_k_coef_temp) %*% f_k_S_temp %*% f_k_coef_temp)
        }
        #if(! (object$type %in% c("FlexGAM1","FlexGAM2"))) warning("No penalty for f_k implemented")
        
        gam_penalty <- sum(gam_penalty_vec)
        # if(norm_beta) {
        #     gam_penalty <- gam_penalty * norm_beta_k
        # }
        dev <- dev + (gam_penalty + f_k_penalty) # 0.5*
    }
    
    as.numeric(dev)
}


print.flexgam <- function(x, ...){
    cat("flexgam object of type:   ", x$type ," \n")
    print(x$family)
    cat("Formula:\n")
    if (is.list(x$formula)){
        for (i in 1:length(x$formula)) {
            print(x$formula[[i]])
        }
    } else {
        print(x$formula)
    }
    
}

summary.flexgam <- function(object, ...){
    print.flexgam(object)
    
    cat("\n Coefficients of predictor:\n")
    coef_mat1 <- cbind(object$beta_k,object$se[1:length(object$beta_k)])
    colnames(coef_mat1) <- c("coef","se")
    rownames(coef_mat1) <- names(object$beta_k)
    print(coef_mat1)
    cat("\n Coefficients of response function:\n")
    coef_mat2 <- cbind(object$f_k$coefficients,object$se[(length(object$beta_k)+1):length(object$se)])
    colnames(coef_mat2) <- c("coef","se")
    rownames_Psi <- names(object$se[(length(object$beta_k)+1):length(object$se)])
    rownames_Psi <- gsub(rownames_Psi,pattern="f_k",replacement="Psi",fixed=TRUE)
    rownames(coef_mat2) <- rownames_Psi
    print(coef_mat2)
}

plot.flexgam <- function(x, type=c("response","covariate"), ci = TRUE,  rug = TRUE, ...){
    type <- match.arg(type)
    dots <- list(...)
    if(is.null(x$F_mat)) ci <- FALSE
    
    if(!is.logical(rug)) stop("rug must be logical")
    if(is.na(rug)) stop("rug must be TRUE/FALSE")
    
    if(!is.logical(ci)) stop("ci must be logical")
    if(is.na(ci)) stop("ci must be TRUE/FALSE")
    
    if("ylim" %in% names(dots)) {
        if(is.list(dots$ylim)) stop("ylim must be a vector of length 2")
        if(is.matrix(dots$ylim)) stop("ylim must be a vector of length 2")
        if(length(dots$ylim) != 2) stop("ylim must be a vector of length 2")
    }
    if("xlim" %in% names(dots)) {
        if(is.list(dots$xlim)) stop("xlim must be a vector of length 2")
        if(is.matrix(dots$xlim)) stop("xlim must be a vector of length 2")
        if(length(dots$xlim) != 2) stop("xlim must be a vector of length 2")
    }
    
    if(type == "response") {
        args2 <- dots
        x_temp <- seq(min(x$eta_k),max(x$eta_k), length=100)
        if(x$type == "FlexGAM2" || x$type == "FlexGAM2n") {
            f_temp <- f_temp_l <- predict_my_pcls_gam(f_k=x$f_k, data=data.frame(eta_k=x_temp), type = x$type)
            X_new_temp <-  PredictMat(x$f_k$bgf$sm,data=data.frame(eta_k=x_temp))
        }
        if(x$type == "FlexGAM1" || x$type == "FlexGAM1n") {
            f_temp_l <- predict(x$f_k, newdata=data.frame(eta_k=x_temp),type="link")
            f_temp <- predict(x$f_k, newdata=data.frame(eta_k=x_temp),type="response")
            X_new_temp <-  PredictMat(x$f_k$smooth[[1]],data=data.frame(eta_k=x_temp))
            if(x$type == "FlexGAM1n") {
                X_new_temp <- cbind(1,X_new_temp)
            }
        }
        
        if(ci) {
            ind1 <- which(grepl(x=colnames(as.matrix(x$F_mat)),pattern="Psi.",fixed=TRUE))
            f_ul_temp <- f_temp_l + 2*sqrt(diag(X_new_temp %*% ginv((as.matrix(x$F_mat))[ind1,ind1]) %*% t(X_new_temp)))
            f_ll_temp <- f_temp_l - 2*sqrt(diag(X_new_temp %*% ginv((as.matrix(x$F_mat))[ind1,ind1]) %*% t(X_new_temp)))
            
            if(x$type == "FlexGAM1" || x$type == "FlexGAM1n") {
                f_ul_temp <- x$family$linkinv(f_ul_temp)
                f_ll_temp <- x$family$linkinv(f_ll_temp)
                
            }
        }
        if(!("ylim" %in% names(dots)) || is.null(dots$ylim)) {
            if(x$family$family == "binomial") {
                args2$ylim <- c(0,1)
            }
            if(x$family$family == "poisson") {
                args2$ylim <- c(0,max(x$fitted))
            }
            if(x$family$family == "gaussian") {
                args2$ylim <- c(min(x$fitted),max(x$fitted))
            }
            if(x$family$family == "Gamma") {
                args2$ylim <- c(0,max(x$fitted))
            }
            
        }
        if(!("ylab" %in% names(dots)) || is.null(dots$ylab)) {
            if(x$family$family == "binomial") {
                args2$ylab <- expression(pi)
            }
            if(x$family$family %in% c("poisson","gaussian","Gamma")) {
                args2$ylab <- expression(mu)
            }
        }
        if(!("xlab" %in% names(dots)) || is.null(dots$xlab)) {
            args2$xlab <- expression(eta)
        }
        args2$type <- "l"
        args2$x <- x_temp
        args2$y <- f_temp
        
        
        do.call(plot, args2 )
        
        if(ci) {
            lines(x_temp, f_ul_temp, lty=2)
            lines(x_temp, f_ll_temp, lty=2)
        }
        if(rug) {
            rug(x$eta_k)
        }
    }
    if(type == "covariate") {
        gp <- interpret.gam(x$formula)
        fake.names <- gp$fake.names
        if(!("ylim" %in% names(dots)) || is.null(dots$ylim)) {
            ylim1 <- c(min(x$eta_k),max(x$eta_k))
        } else {
            ylim1 <- dots$ylim 
        }
        
        if(length(gp$smooth) > 0) {
            vec_smooth_type <- rep(NA,times=length(gp$smooth))
            for(ii in 1:length(gp$smooth)) {
                vec_smooth_type[ii] <- attr(gp$smooth[[ii]],"class")
            }
            plot_number <- which(vec_smooth_type == "ps.smooth.spec")
            if(length(plot_number) > 0) {
                for(ii in plot_number) {
                    args2 <- dots
                    args2$ylim <- ylim1
                    resp_s <- gp$smooth[[ii]]$label
                    resp   <- gp$smooth[[ii]]$term
                    
                    data_temp1 <- x$gam_temp$model[[resp]]
                    x_temp <- seq(min(data_temp1),max(data_temp1),length=100)
                    data_temp <- data.frame(x_temp)
                    names(data_temp) <- resp
                    
                    X_temp <- PredictMat(x$gam_temp$smooth[[ii]],data_temp)
                    f_temp <- X_temp %*% x$beta_k[grepl(names(x$beta_k),pattern=resp,fixed=TRUE)]
                    if(ci){
                        ind2 <- grepl(colnames(as.matrix(x$F_mat)),pattern=resp,fixed=TRUE)
                        f_ul_temp <- f_temp + 2*sqrt(diag(X_temp %*% ginv((as.matrix(x$F_mat))[ind2,ind2]) %*% t(X_temp)))
                        f_ll_temp <- f_temp - 2*sqrt(diag(X_temp %*% ginv((as.matrix(x$F_mat))[ind2,ind2]) %*% t(X_temp)))
                    }
                    ylab1 <- paste("s(",resp,")",sep="")
                    if("ylab" %in% names(dots)) {
                        if(length(dots$ylab) == length(gp$smooth) || length(dots$ylab) == 1) {
                            ylab1 <- dots$ylab[ii]
                        } else {
                            stop("length of ylab doesn't fit to number of smooth effects")
                        }
                    }    
                    xlab1 <- resp
                    if("xlab" %in% names(dots)) {
                        if(length(dots$xlab) == length(gp$smooth) || length(dots$xlab) == 1) {
                            xlab1 <- dots$xlab[ii]
                        } else {
                            stop("length of xlab doesn't fit to number of smooth effects")
                        }
                    }    
                    args2$xlab <- xlab1
                    args2$ylab <- ylab1
                    args2$x <- x_temp
                    args2$y <- f_temp
                    args2$type <- "l"
                    
                    do.call(plot, args2 )
                    
                    if(ci) {
                        lines(x_temp,f_ul_temp,lty=2)
                        lines(x_temp,f_ll_temp,lty=2)
                    }
                    if(rug) {
                        rug(data_temp1)
                    }
                    
                }
            }
            if(length(plot_number) < length(gp$smooth)) {
                warning("Only P-splines of type 'ps' are plotted with estimated confidence intervals. Since here some smooth effects of other type are estimated, their effects are plotted without confidence interval based on the standard output of mgcv.  The P-spline plots are printed twice.")
                args3 <- dots
                args3$ylim <- ylim1
                args3$se <- FALSE
                args3$x <- x$gam_temp
                do.call(plot,args3)
            }
        }
    }
}


response <- function (object, ...) {
    UseMethod("response", object)
}

response.flexgam <- function(object, linear.predictor=NULL, ...) {
    if(is.null(linear.predictor)) linear.predictor <- object$eta_k
    newdata <- data.frame(eta_k=linear.predictor)
    if(object$type == "FlexGAM1" || object$type == "FlexGAM1n") {
        res3 <- predict(object$f_k, newdata=newdata, type="response")
    }
    if(object$type == "FlexGAM2" || object$type == "FlexGAM2n") {
        res3 <- (predict_my_pcls_gam(object$f_k, data=newdata, type = object$type))
    }
    attr(res3,"names") <- NULL
    as.numeric(res3)
    
}

