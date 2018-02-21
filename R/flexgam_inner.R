
Psi_all <- function(eta, f_k, type, family, g_inv = FALSE) {
    if(type == "FlexGAM1n" ) {
        if(g_inv) {
            ret <- predict(f_k,data.frame(eta_k=eta),type="response")
        } else {
            ret <- predict(f_k,data.frame(eta_k=eta),type="link")
        }
    }
    if(type == "FlexGAM1") {
        if(g_inv) {
            ret <- scam::predict.scam(f_k,data.frame(eta_k=eta),type="response",newdata.guaranteed=TRUE)
        } else {
            ret <- scam::predict.scam(object=f_k,newdata=data.frame(eta_k=eta),type="link",newdata.guaranteed=TRUE)
        }
    }
    if(type == "FlexGAM2" || type == "FlexGAM2n") {
        ret <- predict_my_pcls_gam(f_k=f_k, data=data.frame(eta_k=eta), type = type)
        if(family$family %in% c("binomial","poisson","Gamma")) ret[ret <= 0] <- 0
        if(family$family == "binomial") ret[ret >= 1] <- 1
    }
    
    ret
}



dPsi_all <- function(eta, f_k, type, family, Psi_k, deriv_type = "numeric", less_warnings = FALSE) {
    if(deriv_type == "numeric") {
        eps <- 1e-5 ## finite difference interval
        x.mesh <- eta + eps ## shift the evaluation mesh
        X1 <- Psi_all(eta=x.mesh, f_k=f_k, type=type, family = family)
        df <- (X1 - Psi_k)/eps ## maps coefficients to (fd approx.) derivatives
        
    } else {
        if(type == "FlexGAM1") { 
            q2 <- f_k$smooth[[1]]$df + 1
            Sig2 <- matrix(0,q2,q2)   
            for (i in 1:q2)  Sig2[i,1:i] <- 1
            
            # scam package uses splineDesign with outer.ok = FALSE which is not nice, 
            # but induces that for the derivative we have to do the same.
            # Therefore values outside the range of the original eta have to be truncated.
            # This means the derivative is set to constant values outside the range.
            # Even if f_k is not linear outside of the range.
            eta[eta < min(f_k$model$eta_k) ] <- min(f_k$model$eta_k) 
            eta[eta > max(f_k$model$eta_k) ] <- max(f_k$model$eta_k)
            
            Xd1 <- splines::splineDesign(f_k$smooth[[1]]$knots, x=eta, ord=f_k$smooth[[1]]$m+2, derivs=1) 
            df  <- Xd1 %*% Sig2 %*% f_k$coefficients.t
            
        } 
        if(type == "FlexGAM1n") {
            Xd1 <- splines::splineDesign(f_k$smooth[[1]]$knots, x=eta, ord=f_k$smooth[[1]]$m[1]+2, outer.ok=TRUE, derivs=1) 
            df  <- Xd1 %*% f_k$coefficients
            
        }
        if(type == "FlexGAM2" || type == "FlexGAM2n") {
            Xd1 <- splines::splineDesign(knots=f_k$bgf$sm$knots, x=eta, ord=4                 , outer.ok=TRUE, derivs=1)
            df  <- Xd1 %*% f_k$coefficients
            
        }
        df <- df[,1]
    }
    df
}





flexgam_inner <- function(formula, data,
                          type        = c("FlexGAM1","FlexGAM2","FlexGAM1n","FlexGAM2n"),
                          family, 
                          split_list = NULL,
                          estimate_se,
                          control = NULL) {
    
    if(is.null(control )) stop("Major error contact maintainer!")
    
    type <- match.arg(type)
    
    
    control_in            <- control
    control_in$type       <- type
    control_in$formula    <- formula
    
    max_iter_in          <- control_in$max_iter_in
    max_iter_out         <- control_in$max_iter_out 
    delta_in             <- control_in$delta_in 
    delta_out            <- control_in$delta_out
    min_mu_k             <- control_in$min_mu_k
    min_Psi_d            <- control_in$min_Psi_d 
    min_increase         <- control_in$min_increase 
    norm_beta            <- FALSE                   # control_in$norm_beta 
    use_penalty          <- TRUE                    # control_in$use_penalty
    deriv_type           <- "theoretic"             # control_in$deriv_type
    step_halving         <- TRUE                    # control_in$step_halving
    delta_halving        <- control_in$delta_halving
    min_iter_halving_in  <- control_in$min_iter_halving_in
    min_iter_halving_out <- control_in$min_iter_halving_out
    opt_function         <- control_in$opt_function
    start_logit          <- control_in$start_logit
    start_grid_lambda    <- control_in$start_grid_lambda
    fix_smooth           <- control_in$fix_smooth
    sm_par_vec           <- control_in$sm_par_vec
    sp_range             <- control_in$sp_range
    max_iter_out_opt   <- control_in$max_iter_out_opt
    max_iter_in_opt    <- control_in$max_iter_in_opt
    reltol_opt         <- control_in$reltol_opt
    less_warnings        <- TRUE #control_in$less_warnings
    quietly              <- control_in$quietly 
    save_step_response   <- control_in$save_step_response 
    initial_model        <- control_in$initial_model
    
    
    if(!fix_smooth) stop("smoothing parameters must be fixed in loops")
    
    fix_lambda <- TRUE
    
    # sp is inner loop, lambda is outer loop
    
    # Define formula of outer model
    formula_gam  <- y ~ s(eta_k, k=20, bs='ps')
    formula_scam <- y ~ s(eta_k, k=20, bs='mpi')
    
    # Define formula of inner model -> formula_new
    # terms2 and formula_save are for checking of similar sp values in inner loop
    formula_orig <- formula
    formula_new <- formula
    
    formula_new <- update.formula(formula_new,y_k~.)
    formula_new <- update.formula(formula_new,.~.-1)
    formula_i   <- formula_orig
    formula_X_gam <- update.formula(formula_i,.~.-1)
    
    
    
    # Get design matrix for inner loop
    gp <- interpret.gam(formula) # interpret the formula 
    
    y <- model.extract(model.frame(gp$fake.formula,data=data), "response")
    
    n <- nrow(data)
    data$y <- y
    
    data_pred_norm <- data[1:100,]
    if(norm_beta) {
        for(j in 1:ncol(data_pred_norm)) {
            if(!is.numeric(data_pred_norm[,j])) stop("wrong data type, norm beta")
            data_pred_norm[,j] <- seq(min(data[,j]), max(data[,j]) , length = 100)
        }
    }
    
    sm_par_vec_short <- sm_par_vec[which(names(sm_par_vec) != "lambda")]
    if(length(sm_par_vec_short) == 0) sm_par_vec_short <- NULL
    
    # Get start values, gam_temp and beta are for inner loop // p for outer loop
    gam_orig <- gam_temp <- try(gam(formula_i,family=family,data=data,sp=sm_par_vec_short),TRUE)
    gam_X_gam <- try(gam(formula_X_gam,family=family,data=data,sp=sm_par_vec_short, drop.intercept=TRUE),TRUE)
    
    if(inherits(gam_orig, "try-error")){
        gam_orig <- gam_temp <- try(gam(formula_i,  family=gaussian(), data=data,sp=sm_par_vec_short),TRUE)
        warning("Initial model crashed! Applied family=gaussian() in the initial model instead.")
        if(inherits(gam_orig, "try-error")){
            stop("Model crashed at initial values!")
        }
    }
    if(inherits(gam_X_gam, "try-error")){
        gam_X_gam <- try(gam(formula_X_gam,  family=gaussian(), data=data,sp=sm_par_vec_short, drop.intercept=TRUE),TRUE)
        warning("Initial model crashed! Applied family=gaussian() in the initial model instead.")
        if(inherits(gam_X_gam, "try-error")){
            stop("Model crashed at initial values!")
        }
    }
        
    if(initial_model == "with_intercept"){
        beta_orig <- gam_orig$coefficients[-1]
    } else {
        beta_orig <- gam_X_gam$coefficients
    }
    X_gam     <- predict(gam_X_gam, type="lpmatrix")
    
        
    
    norm_beta_k <- sqrt(sum(beta_orig^2))
    
    dev_old <- gam_orig$deviance
    sums_all_f1 <- colMeans_f1 <- sums_all_f <- colMeans_f <- NULL
    
    
    if(!norm_beta){
        beta_temp <- beta_orig
        
        # Center eta for indentifiablity
        eta_orig   <- (X_gam %*% beta_orig)[,1]
        mean_eta_k <- mean(eta_orig)
        sd_eta_k   <- sd(eta_orig)
        
        eta_k    <- eta_orig <- (eta_orig - mean(eta_orig) ) / sd(eta_orig) 
    } else {
        beta_temp <- beta_orig
        eta_orig   <- (X_gam %*% beta_orig)[,1]
        
        # Center eta for indentifiablity
        effects_per_f <- predict(gam_temp, type = "terms")
        
        sums_all_f <- sqrt(sum(predict(gam_temp, type = "terms", newdata=data_pred_norm)^2))
        
        new_effects_per_f_temp <- (effects_per_f)/sums_all_f
        colMeans_f <- colMeans(new_effects_per_f_temp)
        new_effects_per_f  <- new_effects_per_f_temp #- colMeans_f
        
        eta_k <- eta_orig <- rowSums(new_effects_per_f)
        mean_eta_k <- mean(eta_k)
        sd_eta_k   <- sd(eta_k)
        
    }
    
    
    beta_k_old <- beta_k <- beta_temp
    beta_k_old[] <- 0
    beta_list <- list()
    
    
    
    # save coefficients of outer loop    
    p_matrix <- NULL
    p_old  <- NULL
    
    # save outer smoothing parameter and check if outer smoothing paramter is pre-specified   
    
    sp_scam <- lambda_old <- sm_par_vec["lambda"]
    
    
    
    # Save eta, fitted values for plotting of link function at each outer iteration
    eta_matrix    <- NULL
    fitted_matrix <- NULL
    Psi_matrix    <- NULL
    if(save_step_response) {
        eta_matrix    <- eta_orig
        fitted_matrix <- gam_temp$fitted
        Psi_matrix    <- gam_temp$linear.predictor
    }
    
    dev_list <- list()
    
    
    
    # Save Extra Criteria
    
    why_conv_out  <- NULL
    why_conv_in   <- NULL
    where_conv_out  <- NULL
    where_conv_in   <- NULL
    
    freq_wm0 <- NULL
    freq_wm1 <- NULL
    freq_Pd0 <- NULL
    
    # break criteria of outer loop   
    iter_out       <- 1
    conv_out       <- FALSE
    conv_out_dev   <- FALSE
    conv_out_first <- FALSE
    
    while(iter_out <= max_iter_out & !conv_out & !conv_out_dev & !conv_out_first){
        if(!quietly) cat("\n  outer:", iter_out, "   ")
        new_l1 <- FALSE
        if(iter_out == 1) {
            dat <- data.frame(eta_k=eta_k,y=y)
        } else {
            dat$eta_k <- eta_k
        }
        
        # Estimate outer link function
        if(type == "FlexGAM2" || type == "FlexGAM2n") {
            f_k1 <- my_pcls_gam(data=dat, formula=formula_gam, lambda=lambda_old, 
                                family = family, 
                                min_increase=min_increase, 
                                min_mu_k = min_mu_k, 
                                type = type)
        }
        if(type == "FlexGAM1n") {
            f_k1 <- gam(formula=formula_gam, family=family, optimizer=c("outer","nlm"), sp=sp_scam) 
        }
        if(type == "FlexGAM1") {
            f_k1 <- scam::scam(formula=formula_scam, family=family, sp=sp_scam, 
                         data=dat, optimizer = "bfgs", mustart = gam_orig$fitted)
        }
        
        if(type == "FlexGAM1" && 
           ((family$family == "binomial" && ((max(f_k1$fitted) < 1e-15 || min(f_k1$fitted) > 1-1e-15) || any(abs(f_k1$linear.predictor) > 1e5 ))) || 
            (family$family == "poisson" && ((max(f_k1$fitted) < 1e-15 || any(abs(f_k1$linear.predictor) > 1e5 )))) || 
            (family$family == "gaussian" && ( any(abs(f_k1$linear.predictor) > 1e10 ))) || 
            (family$family == "Gamma" && ((max(f_k1$fitted) < 1e-15 || any(abs(f_k1$linear.predictor) > 1e5 ))))) 
        ){
            stop("SCAM optimizing lambda crashed => fixed to break")
        }
        
        
        # Save the values for replotting after finishing the whole process
        if(save_step_response) {
            eta_matrix     <- rbind(eta_matrix    , eta_k)
            fitted_matrix  <- rbind(fitted_matrix , f_k1$fitted)
            Psi_matrix     <- rbind(Psi_matrix    , f_k1$linear.predictor)
        }
        
        
        if(!quietly) cat( "\n  inner: ")
        
        # Save coefficients of link functions
        p <- f_k1$coefficients
        if(is.matrix(p) && dim(p)[1] == 20 && dim(p)[2] == 1) {
            p <- p[,1]
        }
        
        if(save_step_response) {
            
            p_matrix <- rbind(p_matrix,p)
            
            
            
            # Initialize saving of beta coefficients per inner loop
            beta_list[[iter_out]] <- 0
            beta_temp <- NULL
        }
        
        # p_old  is for testing if coefficients of i and i-1 are equal
        if(iter_out <= 1) {
            p_old  <- rep(0,times=length(p))
        }
        
        dev_list[[iter_out]] <- 0
        dev_temp <- NULL
        
        
        
        # break criteria of inner loop   
        iter_in <- 1
        conv_in <- FALSE
        
        if(iter_out == 1){
            dev_k <- dev_new <- deviance_in_fit(eta=eta_k,f_k=f_k1,type=type,y=data$y, 
                                                gam_temp=gam_temp, beta_k = beta_k, 
                                                use_penalty=use_penalty, norm_beta=norm_beta, 
                                                norm_beta_k = norm_beta_k)
        } else {
            dev_k <- dev_new <- deviance_in_fit(eta=eta_k,f_k=f_k1,type=type,y=data$y, 
                                                gam_temp=GG_temp, beta_k = beta_k, 
                                                use_penalty=use_penalty, norm_beta=norm_beta, 
                                                norm_beta_k = norm_beta_k)
        }
        
        # Check if estimated coefficients of link function i and i-1 are equal --> conv_out
        if(sqrt(sum((p_old - p)^2))/sqrt(sum((p^2))) <= delta_out || sqrt(sum((p_old - p)^2)) <= delta_out) {
            conv_out <- TRUE
            why_conv_out <- "conv_out"
            if(!quietly) cat("conv_out \n\r")
        } else {
            # Check if deviance improved, if not stop.
            if(step_halving && min_iter_halving_out <= iter_out && 
               ((dev_new-dev_old)/(0.1 + abs(dev_new)) >= -delta_halving)) {
                why_conv_out <- "conv_out_dev"
                conv_out_dev <- TRUE
                conv_out_first <- FALSE
                if(iter_out == 1) {
                    f_k <- f_k1
                    conv_out_first <- TRUE
                    conv_out_dev   <- FALSE
                    why_conv_out   <- "conv_out_first"
                    if(!quietly) cat("conv_out_first \n\r")
                } else {
                    if(!quietly) cat("conv_out_dev \n\r")
                }
            } else {
                conv_out_dev <- FALSE
                conv_out_first <- FALSE
                why_conv_out <- "no"
                #dev_old <- dev_new
                f_k <- f_k1
                p_old <- p
                
            }
        } # End check for conv_out
        
        while(iter_in <= max_iter_in & !conv_in & !conv_out & !conv_out_dev & !conv_out_first) {
            new_l1 <- FALSE
            
            if(!quietly) {
                cat(iter_in, " ")
                if(iter_in %% 20 == 0)   cat("\n\r")
            }
            
            dev_old <- dev_k
            
            
            # Estimate mu_k = pi_k
            # Estimate derivative of estimated link function at the current eta
            Psi_k  <-   Psi_all(eta=eta_k,f_k=f_k,type=type, family = family)
            Psi_d  <-  dPsi_all(eta=eta_k,f_k=f_k,type=type, family = family, 
                                Psi_k = Psi_k, deriv_type = deriv_type, less_warnings = less_warnings)
            
            if(type == "FlexGAM1" || type == "FlexGAM1n"){
                mu_k <- family$linkinv(Psi_k)
            } else {
                mu_k <- Psi_k
            }
            
            # Very small values of mu and dPsi cause crashs, prevent them
            if(family$family %in% c("binomial","poisson","Gamma")) {
                if(any(mu_k < min_mu_k)){
                    wm0 <- which(mu_k < min_mu_k)
                    mu_k[wm0] <- min_mu_k
                    freq_wm0 <- c(freq_wm0,length(wm0)/length(mu_k))
                } else {
                    freq_wm0 <- c(freq_wm0,0)
                }
            }
            if(family$family == "binomial") {
                if(any(mu_k > 1 - min_mu_k)){
                    wm1 <- which(mu_k > 1 - min_mu_k)
                    mu_k[wm1] <- 1-min_mu_k
                    freq_wm1 <- c(freq_wm1,length(wm1)/length(mu_k))
                } else {
                    freq_wm1 <- c(freq_wm1,0)
                }
            }
            if(any(Psi_d < min_Psi_d)){
                Pd0 <- which(Psi_d < min_Psi_d)
                Psi_d[Pd0] <- min_Psi_d
                freq_Pd0 <- c(freq_Pd0,length(Pd0)/length(Psi_d))
            } else {
                freq_Pd0 <- c(freq_Pd0,0)
            }
            
            # Build y_k and w_k
            if(family$family == "binomial") {
                if(type == "FlexGAM1" || type == "FlexGAM1n") {
                    y_k <- eta_k + (y - mu_k) / (mu_k * (1 - mu_k) * Psi_d) 
                    w_k <- mu_k * (1 - mu_k) * (Psi_d)^2
                } else {
                    y_k <- (eta_k + (y - mu_k) / (Psi_d) )
                    w_k <- (Psi_d^2/(mu_k*(1-mu_k)))
                }
            } 
            if(family$family == "poisson") {
                if(type == "FlexGAM1" || type == "FlexGAM1n") {
                    y_k <- eta_k + (y - mu_k) / (mu_k * Psi_d) 
                    w_k <- mu_k * (Psi_d)^2
                } else {
                    y_k <- (eta_k + (y - mu_k) / (Psi_d) )
                    w_k <- ((Psi_d)^2/mu_k)
                }
            }
            if(family$family == "gaussian") {
                if(type == "FlexGAM1" || type == "FlexGAM1n") {
                    y_k <- (eta_k + (y - mu_k) / (Psi_d) )
                    #phi_k <- sum((y-mu_k)^2)/(n - sum(gam_temp$edf))
                    #phi_k <- 1
                    w_k <- ((Psi_d)^2)
                } else {
                    y_k <- (eta_k + (y - mu_k) / (Psi_d) )
                    #phi_k <- sum((y-mu_k)^2)/(n - sum(gam_temp$edf))
                    #phi_k <- 1
                    w_k <- ((Psi_d)^2)
                }
            }
            if(family$family == "Gamma") {
                if(type == "FlexGAM1" || type == "FlexGAM1n") {
                    y_k <- eta_k + ((y - mu_k) / (mu_k*Psi_d ))
                    w_k <- ((Psi_d)^2)
                } else {
                    y_k <- (eta_k + (y - mu_k) / (Psi_d) )
                    w_k <- ((Psi_d)^2/mu_k^2)
                }
            }
            if(any(Psi_d < 1e-8)) {
                w_k[Psi_d < 1e-8] <- 0
                warning(paste0(type, " Pd0_2",sep=" "))
            }
            
            if(iter_out == 1 && iter_in == 1) {
                data_new <- data.frame(data,y_k=y_k,w_k=w_k)
            } else {
                data_new$y_k <- y_k
                data_new$w_k <- w_k
            }
            
            # Estimate the effects of the true covariates
            #gam_temp  <- gam(formula_new, data=data_new, weights=w_k, optimizer=c("outer","nlm"), fit=TRUE)
            if(iter_in == 1 && iter_out == 1) {
                GG_temp  <- gam(formula_new, data=data_new, weights=w_k, fit=FALSE, sp=sm_par_vec_short, drop.intercept=TRUE)
                
            } else {
                GG_temp$y <- y_k
                GG_temp$w <- w_k
                GG_temp$mf$y_k <- y_k
                GG_temp$mf[["(weights)"]] <- w_k
            }
            
            gam_temp1  <- gam(G=GG_temp, fit=TRUE)
            beta_k1 <- gam_temp1$coefficients
            norm_beta_k1 <- sqrt(sum(beta_k1^2))
            
            if(!norm_beta) {
                eta_k1 <- (X_gam %*% beta_k1)[,1]
                mean_eta_k1 <- mean(eta_k1)
                sd_eta_k1   <- sd(eta_k1)
                eta_k1 <- ( eta_k1 - mean(eta_k1) ) / sd(eta_k1)
            } else {
                
                # Center eta for indentifiablity
                effects_per_f <- predict(gam_temp1, type = "terms")
                
                sums_all_f1 <- sqrt(sum(predict(gam_temp1, type = "terms", newdata=data_pred_norm)^2))
                
                new_effects_per_f_temp <- (effects_per_f)/sums_all_f1
                colMeans_f1 <- colMeans(new_effects_per_f_temp)
                new_effects_per_f  <- new_effects_per_f_temp #- colMeans_f1
                
                eta_k1 <- rowSums(new_effects_per_f)
                
                mean_eta_k1 <- mean(eta_k1)
                sd_eta_k1   <- sd(eta_k1)
                
            }    
            
            dev_new <- deviance_in_fit(eta=eta_k1, f_k=f_k, type=type, y=data$y, 
                                       gam_temp=GG_temp, beta_k = beta_k1,
                                       use_penalty=use_penalty, norm_beta=norm_beta, 
                                       norm_beta_k = norm_beta_k1)
            
            conv_in_halv <- FALSE
            if(step_halving && (iter_out >= 2 || iter_in >= 2) && min_iter_halving_in <= iter_in  ){   
                iter_in2 <- 1
                new_l1 <- FALSE
                while((dev_new - dev_old)/(0.1 + abs(dev_new)) > -delta_halving && iter_in2 <= 20) {
                    new_l1 <- TRUE
                    if(!quietly) cat(paste("s", iter_in2, " ",sep=""))
                    beta_k1 <- (beta_k1 + beta_k)/2	
                    
                    if(!norm_beta) {
                        gam_temp1$coefficients <- beta_k1
                        eta_k1   <- (X_gam %*% beta_k1)[,1]
                        mean_eta_k1 <- mean(eta_k1)
                        sd_eta_k1   <- sd(eta_k1)
                        eta_k1 <- ( eta_k1 - mean(eta_k1) ) / sd(eta_k1)
                    } else {
                        
                        gam_temp1$coefficients <- beta_k1
                        
                        # Center eta for indentifiablity
                        effects_per_f <- predict(gam_temp1, type = "terms")
                        
                        sums_all_f1 <- sqrt(sum(predict(gam_temp1, type = "terms", newdata=data_pred_norm)^2))
                        
                        new_effects_per_f_temp <- (effects_per_f)/sums_all_f1
                        colMeans_f1 <- colMeans(new_effects_per_f_temp)
                        new_effects_per_f  <- new_effects_per_f_temp # - colMeans_f1
                        
                        eta_k1 <- rowSums(new_effects_per_f)
                        
                        mean_eta_k1 <- mean(eta_k1)
                        sd_eta_k1   <- sd(eta_k1)
                    }    
                    
                    
                    
                    
                    dev_new <- deviance_in_fit(eta=eta_k1, f_k=f_k, type=type, y=data$y, 
                                               gam_temp=GG_temp, beta_k = beta_k1, 
                                               use_penalty=use_penalty, norm_beta=norm_beta, 
                                               norm_beta_k = norm_beta_k1)
                    
                    iter_in2 <- iter_in2 + 1
                }
                if(iter_in2 > 20) {
                    beta_k1 <- beta_k
                    eta_k1 <- eta_k
                    mean_eta_k1 <- mean_eta_k
                    sd_eta_k1 <- sd_eta_k
                    colMeans_f1 <- colMeans_f
                    sums_all_f1 <- sums_all_f
                    conv_in_halv <- TRUE
                    dev_k <- dev_old
                    norm_beta_k1 <- norm_beta_k
                    gam_temp1 <- gam_temp
                } else {
                    dev_k <- dev_new
                    dev_old <- dev_new
                }
            } else {
                dev_k <- dev_new
            }
            if(!quietly && new_l1) cat("\n         ")
            
            if(length(beta_k) != length(beta_k1)) stop("length(beta_k) != length(beta_k1)")
            
            # Check if estimated coefficients of i an i-1 are equal --> convergence (conv_in)
            if(sqrt(sum((beta_k - beta_k1)^2))/sqrt(sum((beta_k1^2))) <= delta_in || sqrt(sum((beta_k - beta_k1)^2)) <= delta_in){ 
                conv_in <- TRUE
                if(!conv_in_halv) {
                    why_conv_in <- c(why_conv_in,"conv_in")
                    if(!quietly) {cat(" conv_in ")}
                } else {
                    why_conv_in <- c(why_conv_in,"conv_in_halv")
                    if(!quietly) {cat(" conv_in_halv ")}
                }
                where_conv_in <- c(where_conv_in,iter_in)
            } else {
                conv_in <- FALSE
                beta_k_old <- beta_k
                beta_k     <- beta_k1
                norm_beta_k <- norm_beta_k1
                eta_k      <- eta_k1
                mean_eta_k <- mean_eta_k1
                sd_eta_k   <- sd_eta_k1
                colMeans_f <- colMeans_f1
                sums_all_f <- sums_all_f1
                gam_temp <- gam_temp1
                
            } # End check for conv_in
            
            
            # Save estimated coefficients of inner model
            if(save_step_response) {
                beta_temp <- rbind(beta_temp,beta_k)
                dev_temp <- c(dev_temp,dev_k)
            }
            
            # Save smoothing parameter of inner model
            
            iter_in <- iter_in + 1
            if(iter_in == max_iter_in + 1 && !conv_in) {
                why_conv_in <- c(why_conv_in,"no")
                where_conv_in <- c(where_conv_in,iter_in)
            }
        } # End inner loop
        
        # Save parameters of inner loop
        if(save_step_response) {
            beta_list[[iter_out]] <- beta_temp
            dev_list[[iter_out]] <- dev_temp
        }
        
        
        iter_out <- iter_out + 1
        # if(iter_out == max_iter_out) {
        #     if(!less_warnings) warning(paste0(type, " max_iter_out reached",sep=" "))
        #     if(iter_in == max_iter_in) warning(paste0(type, "max_iter_in reached",sep=" "))
        # }
        
        if(!quietly)    cat("\n")
        
    } # End outer loop
    
    
    where_conv_out <- iter_out
    
    # Estimate fitted values
    Psi_k  <-   Psi_all(eta=eta_k, f_k=f_k, type=type, family = family)
    Psi_d  <-  dPsi_all(eta=eta_k, f_k=f_k, type=type, family = family, 
                        Psi_k = Psi_k, deriv_type = deriv_type, less_warnings = less_warnings)
    
    if(type == "FlexGAM1" || type == "FlexGAM1n") { 
        fitted <- family$linkinv(Psi_k)
    } else {
        fitted <- Psi_k
    }
    
    if(length(where_conv_in) == 1 && where_conv_in == 1) conv_out_first <- TRUE
    
    # Save fitted values of inner model for plotting
    if(conv_out_first) { ## iter_out was updated to 2, even, if conv_out_dev is TRUE
        if(save_step_response) {
            beta_list   [[1]] <- as.matrix(beta_k,nrow=1,ncol=length(beta_k))
            dev_list    [[1]] <- gam_orig$deviance
        }
        
        why_conv_in <- "no_inner_loop"
        where_conv_in <- NA
        freq_wm0 <- 0
        freq_wm1 <- 0
        freq_Pd0 <- 0
        
        
    }
    
    
    if(save_step_response) {
        deviance <- unlist(dev_list)[length(unlist(dev_list))]
    } else {
        deviance <- deviance_in_fit(eta=eta_k, f_k=f_k, type=type, y=data$y, 
                                    gam_temp=GG_temp, beta_k = beta_k, 
                                    use_penalty=use_penalty, norm_beta=norm_beta, 
                                    norm_beta_k = norm_beta_k)
    }
    
    sp_ret <- gam_temp$sp
    if(is.null(sp_ret) || length(gam_temp$sp) == 0) {
        if(length(gam_temp$smooth) > 0) {
            sp_ret <- NULL
            for(ispr in 1:length(gam_temp$smooth)) {
                sp_ret <- c(sp_ret,gam_temp$smooth[[ispr]]$sp)
            }
        }
    }
    
    phi <- 1
    
    if(family$family == "gaussian") {
        phi <- sum((y-fitted)^2)/(length(y) - sum(gam_temp$edf))
    } 
    if(family$family == "Gamma") {
        phi <- sum(((y-fitted)^2)/((length(y) - sum(gam_temp$edf))*(fitted^2)))
    } 
    
    se1 <- NULL
    if(estimate_se){
        se1 <- build_se(type=type, gam_temp=gam_temp, f_k=f_k, eta_k=eta_k,
                        Psi_k=Psi_k, Psi_d=Psi_d, fitted=fitted, family=family, phi=phi) 
    }
    
    coefficients <- c(beta_k,f_k$coefficients)
    names(coefficients) <- c(names(beta_k), paste("Psi.",1:length(f_k$coefficients),sep=""))
    
    ret1 <- list("gam_orig"       = gam_orig, 
                 "response"       = y,
                 "fitted"         = fitted, 
                 "linear.predictor" = eta_k, 
                 "coefficients"   = coefficients,
                 "formula"        = formula,
                 "type"           = type, 
                 "family"         = family,
                 "control"        = control_in, 
                 "split_list"     = split_list,
                 "beta_k"         = beta_k,
                 "se"             = se1$se,
                 "F_mat"          = se1$F_mat,
                 "norm_beta_k"    = norm_beta_k, 
                 "eta_k"          = eta_k, 
                 "f_k"            = f_k, 
                 "Psi_k"          = Psi_k,
                 "Psi_d"          = Psi_d,
                 "mean_eta_k"     = mean_eta_k, 
                 "sd_eta_k"       = sd_eta_k,
                 #"colMeans_f"     = colMeans_f,
                 #"sums_all_f"     = sums_all_f,
                 "gam_temp"       = gam_temp,
                 "sm_par_vec"     = sm_par_vec,
                 "fix_smooth"     = fix_smooth, 
                 "formula_new"    = formula_new,
                 "deviance"       = deviance,
                 "data_pred_norm" = data_pred_norm,
                 "details" = list(
                     "why_conv_in"    = why_conv_in   , "where_conv_in"    = where_conv_in, 
                     "why_conv_out"   = why_conv_out  , "where_conv_out"   = where_conv_out,
                     "freq_wm0"       = freq_wm0      , "freq_wm1"         = freq_wm1,
                     "freq_Pd0"       = freq_Pd0,
                     "eta_matrix"     = eta_matrix, 
                     "fitted_matrix"  = fitted_matrix, 
                     "Psi_matrix"     = Psi_matrix,
                     "iter_out"       = iter_out, 
                     "iter_in"        = iter_in,
                     "conv_in"        = conv_in, 
                     "conv_in_halv"   = conv_in_halv, 
                     "conv_out"       = conv_out,  
                     "conv_out_dev"   = conv_out_dev, 
                     "conv_out_first" = conv_out_first,
                     "p_matrix"       = p_matrix,
                     "beta_list"      = beta_list, 
                     "lambda"         = sp_scam, 
                     "sp"             = sp_ret,
                     "dev_list"       = dev_list)
    )
    class(ret1) = c("flexgam")
    ret1
}



deviance_in_fit <- function(eta, f_k, type, y, gam_temp, beta_k, use_penalty, norm_beta, norm_beta_k) {
    if(type == "FlexGAM2" || type == "FlexGAM2n") {
        fitted <- Psi_all(eta=eta,f_k=f_k,type=type, family = f_k$family)
    }
    if(type == "FlexGAM1" || type == "FlexGAM1n") {
        fitted <- Psi_all(eta=eta,f_k=f_k,type=type, family = f_k$family, g_inv = TRUE)
    }
    
    if(f_k$family$family == "binomial" && any(fitted >= 1-1e-16)) fitted[fitted >= 1-1e-16] <- 1-1e-16 
    # Since predict.gam uses linear extrapolation it might occure 
    if(f_k$family$family %in% c("binomial","poisson","Gamma") && any(fitted <    1e-16)) fitted[fitted <    1e-16] <-   1e-16 
    # that the predicted values are outside the range of [0,1].
    
    if(f_k$family$family == "binomial")  dev_k <- -2*(sum(log(1-fitted[y==0])) + sum(log(fitted[y==1])))
    if(f_k$family$family == "poisson" )  dev_k <- 2*(sum(y[y>0] * log(y[y>0] / fitted[y>0])) + sum( - (y - fitted)))
    if(f_k$family$family == "gaussian" ) dev_k <- sum((y - fitted)^2)
    if(f_k$family$family == "Gamma" )    dev_k <- 2*sum((y - fitted)/fitted - log(y/fitted))
    
    if(use_penalty) {
        gam_penalty_vec <- 0
        if(length(gam_temp$smooth) > 0) {
            gam_penalty_vec <- NULL
            for(jj in 1: length(gam_temp$smooth)){
                gam_S_temp  <- gam_temp$smooth[[jj]]$S
                gam_sp_temp <- gam_temp$smooth[[jj]]$sp
                gam_S_temp1 <- gam_S_temp[[1]]
                gam_S_temp1[] <- 0
                for(jj_1 in 1:length(gam_S_temp)) {
                    gam_S_temp1 <- gam_S_temp1 + gam_S_temp[[jj_1]] * gam_sp_temp[jj_1]
                }
                if(any(gam_sp_temp == -1)) stop("Major error sp should be fixed!")
                gam_coef_temp <- beta_k[grepl(x=names(beta_k),pattern=gam_temp$smooth[[jj]]$label,fixed=TRUE)]
                #gam_coef_temp <- beta_k[1:19 + (19*(jj-1))]
                gam_penalty_vec <- c(gam_penalty_vec, (t(gam_coef_temp) %*% gam_S_temp1 %*% gam_coef_temp))
            }
        }
        f_k_penalty <- 0
        if(type == "FlexGAM1"  || type == "FlexGAM1n"){
            f_k_S_temp  <- f_k$smooth[[1]]$S[[1]]
            f_k_sp_temp <- f_k$smooth[[1]]$sp
            if(f_k_sp_temp == -1) f_k_sp_temp <- f_k$sp[names(f_k_sp_temp)]
            f_k_coef_temp <- f_k$coefficients[grepl(x=names(f_k$coefficients),pattern=names(f_k_sp_temp),fixed=TRUE)]
            f_k_penalty <- f_k_sp_temp * (t(f_k_coef_temp) %*% f_k_S_temp %*% f_k_coef_temp)
        }
        if(type == "FlexGAM2" || type == "FlexGAM2n") {
            f_k_S_temp  <- f_k$bgf$sm$S[[1]]
            f_k_sp_temp <- f_k$sp
            f_k_coef_temp <- f_k$coefficients
            f_k_penalty <- f_k_sp_temp * (t(f_k_coef_temp) %*% f_k_S_temp %*% f_k_coef_temp)
        }
        #if(! (type %in% c("FlexGAM1","FlexGAM2"))) warning("No penalty for f_k implemented")
        
        gam_penalty <- sum(gam_penalty_vec)
        # if(norm_beta) {
        #     gam_penalty <- gam_penalty * norm_beta_k
        # }
        dev_k <- dev_k + (gam_penalty + f_k_penalty) # 0.5*
    }
    as.numeric(dev_k)
}



chunk3 <- function(vec,ng) {
    nv <- length(vec)
    length1 <- nv %/% ng
    length_sample <- rep(length1,times=ng)
    if(nv %% ng > 0) {
        length_sample[1:(nv %% ng)] <- length_sample[1:(nv %% ng)] +1
    }
    
    res_l <- list()
    vec1 <- vec
    for(ichunk in 1:ng) {
        res_l[[ichunk]] <- sample(vec1, size = length_sample[ichunk])
        vec1 <- vec1[ !(vec1 %in% res_l[[ichunk]])]
    }
    res_l
}


AUC_function_short <- function(fitted, y, length_roc = 100) {
    AUC1 <- NA
    if(all(y %in% c(0,1))) {
        TruePositiveRate1 <- FalsePositiveRate1 <- rep(NA, times=length_roc)
        roc_grid <- seq(0,1,length=length_roc)
        for(roc in 1:length_roc) {
            theta <- roc_grid[roc]
            TruePositiveRate1[roc]  <- sum(y == 1 & fitted >= theta) / sum(y==1)
            FalsePositiveRate1[roc] <- sum(y == 0 & fitted >= theta) / sum(y==0)
        }
        diffFPR1 <- -diff(FalsePositiveRate1)
        AUC1a <- sum(TruePositiveRate1[-1]*diffFPR1)
        AUC1b <- sum(TruePositiveRate1[-length(TruePositiveRate1)]*diffFPR1)
        AUC1 <- mean(c(AUC1a,AUC1b))
    } 
    AUC1
}





get_penalty <- function(type,gam_temp,f_k) {
    gam_penalty_vec <- 0
    gam_penalty <- 0
    if(length(gam_temp$smooth) > 0) {
        gam_penalty_vec <- NULL
        for(jj in 1: length(gam_temp$smooth)){
            gam_S_temp  <- gam_temp$smooth[[jj]]$S
            gam_sp_temp <- gam_temp$smooth[[jj]]$sp
            gam_S_temp1 <- gam_S_temp[[1]]
            gam_S_temp1[] <- 0
            for(jj_1 in 1:length(gam_S_temp)) {
                gam_S_temp1 <- gam_S_temp1 + gam_S_temp[[jj_1]] * gam_sp_temp[jj_1]
            }
            if(any(gam_sp_temp == -1)) stop("Major error sp should be fixed!")
            
            gam_penalty <- bdiag(gam_penalty,gam_S_temp1)
        }
        gam_penalty <- gam_penalty[-1,-1]
        if(gam_temp$nsdf > 0){
            for(i in 1:gam_temp$nsdf) {
                gam_penalty <- cbind(0,rbind(0,gam_penalty))
            }
        }
    }
    f_k_penalty <- 0
    if(type == "FlexGAM1" || type == "FlexGAM1n"){
        f_k_S_temp  <- f_k$smooth[[1]]$S[[1]]
        f_k_sp_temp <- f_k$smooth[[1]]$sp
        if(f_k_sp_temp == -1) f_k_sp_temp <- f_k$sp[names(f_k_sp_temp)]
        f_k_penalty <- f_k_sp_temp *  f_k_S_temp 
        f_k_penalty <- cbind(0,rbind(0,f_k_penalty))
    }
    if(type == "FlexGAM2" || type == "FlexGAM2n") {
        f_k_S_temp  <- f_k$bgf$sm$S[[1]]
        f_k_sp_temp <- f_k$sp
        f_k_penalty <- f_k_sp_temp * f_k_S_temp 
    }
    list("gam_penalty"=gam_penalty,"f_k_penalty"=f_k_penalty)
}

build_se <- function(type, gam_temp, f_k, eta_k, Psi_k, Psi_d, fitted, family, phi=NULL) {
    if(family$family %in% c("binomial","poisson","gaussian","Gamma")) {
        if(family$family == "binomial") {
            if(any(fitted >= 1-1e-16)) fitted[fitted >= 1-1e-16] <- 1-1e-16 
            # Since predict.gam uses linear extrapolation it might occure 
             # that the predicted values are outside the range of [0,1].
        }
        if(family$family %in% c("binomial","poisson","Gamma")) {
            if(any(fitted <    1e-16)) fitted[fitted <    1e-16] <-   1e-16 
        }
        
        X_gam <- predict(gam_temp, type="lpmatrix")
        if(type=="FlexGAM1" || type == "FlexGAM1n") {
            B_gam <- predict(f_k,newdata=data.frame(eta_k=eta_k), type="lpmatrix")
        }
        if(type=="FlexGAM2" || type == "FlexGAM2n") {
            B_gam <- f_k$bgf$M$X
        }
        Fbb <- matrix(0,nrow=ncol(X_gam),ncol=ncol(X_gam))
        Fnn <- matrix(0,nrow=ncol(B_gam),ncol=ncol(B_gam))
        Fbn <- matrix(0,nrow=ncol(B_gam),ncol=ncol(X_gam))
        hh_gam <- NULL
        if(family$family == "binomial") {
            if(type=="FlexGAM1" || type == "FlexGAM1n") {
                hh_gam <- fitted*(1-fitted)
            }
            if(type=="FlexGAM2" || type == "FlexGAM2n") {
                hh_gam <- 1/(fitted*(1-fitted))
            }
        }
        if(family$family == "poisson") {
            if(type=="FlexGAM1" || type == "FlexGAM1n") {
                hh_gam <- fitted
            }
            if(type=="FlexGAM2" || type == "FlexGAM2n") {
                hh_gam <- 1/fitted
            }
        }
        if(family$family == "gaussian") {
            if(type=="FlexGAM1" || type == "FlexGAM1n") {
                hh_gam <- rep(1/phi, times = nrow(X_gam))
            }
            if(type=="FlexGAM2" || type == "FlexGAM2n") {
                hh_gam <- rep(1/phi, times = nrow(X_gam))
            }
        }
        if(family$family == "Gamma") {
            if(type=="FlexGAM1" || type == "FlexGAM1n") {
                hh_gam <- rep(1/phi, times = nrow(X_gam))
            }
            if(type=="FlexGAM2" || type == "FlexGAM2n") {
                hh_gam <- 1/(fitted^2)/phi
            }
        }
        for(ii in 1:nrow(X_gam)) {
            x1 <- X_gam[ii,]
            B1 <- B_gam[ii,]
            Fbb <- Fbb + x1%*%t(x1)*hh_gam[ii]*(Psi_d[ii])^2 
            Fnn <- Fnn + B1%*%t(B1)*hh_gam[ii] 
            Fbn <- Fbn + B1%*%t(x1)*hh_gam[ii]*(Psi_d[ii]) 
        }
        Fbb <- Fbb + get_penalty(type, gam_temp, f_k)$gam_penalty/phi
        Fnn <- Fnn + get_penalty(type, gam_temp, f_k)$f_k_penalty/phi
        
        F_mat <- rbind(cbind(Fbb,t(Fbn)),cbind(Fbn,Fnn))
        diag_se1 <- diag(ginv(as.matrix(F_mat)))
        if(any(diag_se1 < 0 & diag_se1 > -1e-8)) {
            warning("Penalty too high! Variance of coefficients is numerical 0.")
            diag_se1[diag_se1 < 0 & diag_se1 > -1e-8] <- 0
        }
        se <- sqrt(diag_se1)
        names(se) <- c(names(gam_temp$coefficients),paste("Psi.",1:length(f_k$coefficients),sep=""))
        colnames(F_mat) <- rownames(F_mat) <- c(names(gam_temp$coefficients),paste("Psi.",1:length(f_k$coefficients),sep=""))
    } else {
        se <- NULL
        F_mat <- NULL
    }
    list("se"=se, "F_mat" = F_mat)
}
