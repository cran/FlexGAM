flexgam <- function(formula, 
                    data, 
                    type        = c("FlexGAM1","FlexGAM2","FlexGAM1n","FlexGAM2n"),
                    family,
                    control = NULL
){
    
    type <- match.arg(type)
    cl <- match.call()
    control_input <- control
    if(attr(terms.formula(formula),"intercept") == 0) {
        stop("formula is wrong. Intercept will be excluded automatically!")
    }
    
    if(class(data) != "data.frame") {
        stop("data must be a data.frame")
    }
     
    if(!(family$family %in% c("binomial","poisson","gaussian","Gamma") )) stop("Only binomial, poisson and gaussian as families implemented!")
    if(type == "FlexGAM2n" && family$family == "gaussian") stop("FlexGAM2n not possible with gaussian data. Use FlexGAM1n or FlexGAM2 instead.")
    
    if(length(colnames(data)) > 0 && any(colnames(data) == "y") && as.character(formula)[2] != "y" ) {
        stop("If y is a colname of data it has to be the response")
    }
    
    gp <- interpret.gam(formula) # interpret the formula 
    
    if(!gp$response %in% names(data)) {
        stop("Response missing in data")
    }
    
    if(length(unique(colnames(data))) != length((colnames(data)))) {
        stop("Duplicate colnames exist")
    }
    
    mf <- model.frame(gp$fake.formula,data=data, na.action=na.pass)
    y <- model.extract(mf, "response")
    
    if(!(all(gp$fake.names %in% colnames(data)))) stop("Some variables are not included in the data.")
    
    if(!is.numeric(y)) stop("Response has to be numeric")
    
    if(any(is.na(mf))) stop("Data must not contain NA!")
    if(dim(mf)[2] < 3) stop("There must be at least 2 covariates to fit a FlexGAM model!")
    
    included_var <- c(gp$response,gp$fake.names)
    data <- data[,included_var]
    
    control <- match_flexgam_control(control = control, formula = formula, data=data)
    
    start_logit       <- control$initial_sm
    start_grid_lambda <- control$initial_sm
    fix_smooth        <- control$fix_smooth
    
    if(control$min_iter_halving_out < 2) stop("min_iter_halving_out must be at least 2.")
    
    if(family$family == "binomial" && !all(y %in% c(0,1))) stop("Response must be either 0 or 1")
    
    if(family$family == "Gamma"    && family$link != "log")      stop("Gamma    is only implemented for log link!")
    if(family$family == "poisson"  && family$link != "log")      stop("Poisson  is only implemented for log link!")
    if(family$family == "gaussian" && family$link != "identity") stop("Gaussian is only implemented for identity link!")
    if(family$family == "binomial" && family$link != "logit")    stop("Binomial is only implemented for logit link!")
    
    
    # if(control$deriv_type != "numeric" && control$type != "FlexGAM2") {
    #     warning("deriv_type must be theoretic, if type is != FlexGAM2! Therefore deriv_type is changed to numeric")
    #     control$deriv_type <- "numeric"
    # }
    
    #if(!exists("split_list") || is.null(split_list)) {
    split_list <- chunk3(vec=1:nrow(data),ng = 5)
    #}
    
    if(!fix_smooth){
        if(start_logit) {
            start_model <- gam(formula, data=data, family=family)
            temp_sp <- start_model$sp 
            min_logit_sp <- max(c(1e-6,control$sp_range[1]*10 ))
            max_logit_sp <- min(c(1e+6,control$sp_range[2]/10))
            if(any(temp_sp < min_logit_sp)) temp_sp[temp_sp < min_logit_sp] <- min_logit_sp
            if(any(temp_sp > max_logit_sp)) temp_sp[temp_sp > max_logit_sp] <- max_logit_sp
            control$sm_par_vec <- c(control$sm_par_vec["lambda"],temp_sp)
        }
        
        if(start_grid_lambda){
            min_grid_lambda <- -2
            max_grid_lambda <- 4
            if(ceiling(log10(control$sp_range[1])) > min_grid_lambda) min_grid_lambda <- ceiling(log10(control$sp_range[1]))
            if(floor  (log10(control$sp_range[2])) < max_grid_lambda) max_grid_lambda <- floor  (log10(control$sp_range[2]))
            grid_lambda <- 10^seq(min_grid_lambda,max_grid_lambda,by=1)
            
            control_grid <- control
            control_grid$fix_smooth <- TRUE
            control_grid$quietly    <- TRUE
            
            control_grid$save_step_response <- FALSE
            
            dev_grid <- rep(NA, times=length(grid_lambda))
            for(tt in 1:length(grid_lambda)) {
                control_grid$sm_par_vec["lambda"] <- grid_lambda[tt]
                dev_vec_grid <- rep(NA,times=5)
                for(ig1 in 1:5) {
                    ret_grid <- try(flexgam_inner(formula = formula, 
                                                  data        = data[-split_list[[ig1]],], 
                                                  type        = type,
                                                  family      = family, 
                                                  split_list  = split_list,
                                                  estimate_se = FALSE,
                                                  control = control_grid),TRUE)
                    if(!inherits(ret_grid, "try-error")){
                        dev_vec_grid[ig1] <- deviance(object=ret_grid,newdata=data[split_list[[ig1]],])
                    } else {
                        dev_vec_grid[ig1] <- NA
                    }
                }
                if(!all(is.na(dev_vec_grid))) {
                    dev_grid[tt] <- mean(dev_vec_grid,na.rm=TRUE)
                } 
            }
            if(!control$quietly){
                cat(" Results of initial crossvalidation for lambda:\n")
                temp_grid_res <- rbind(grid_lambda,dev_grid)
                rownames(temp_grid_res) <- c("grid","CV")
                print(temp_grid_res,colnames=FALSE)
            }
            if(!all(is.na(dev_grid))) {
                control$sm_par_vec["lambda"] <- grid_lambda[which.min(dev_grid)]
            } else {
                warning("start grid lambda resulted in NA only!")
            }
        }
        if(!control$quietly){
            cat("\n","Initial values for optimizing all smoothing parameters jointly:\n")
            print(control$sm_par_vec)
            cat("\n","Optimizing trace in log-scale: \n")
        }
        ret1 <- flexgam_optimize(formula_1a  = formula, 
                                 data_1a     = data, 
                                 type        = type,
                                 family      = family, 
                                 split_list  = split_list,
                                 control = control )
        
    } else {
        control_fixed <- control
        control_fixed$fix_smooth <- TRUE
        if(!control_fixed$quietly){
            cat(" Smoothing parameters:\n")
            print(control_fixed$sm_par_vec)
        }
        ret1 <- flexgam_inner(formula     = formula, 
                              data        = data, 
                              type        = type,
                              family      = family, 
                              split_list  = split_list,  
                              estimate_se = TRUE,
                              control = control_fixed)
    }
    ret1$call <- cl
    ret1$control_input <- control_input
    ret1
}





match_flexgam_control <- function(control = NULL, formula = formula, data=NULL) {
    if(is.null(control)) {
        control_new <- list("max_iter_in"        = 100,
                            "max_iter_out"       = 100, 
                            "delta_in"           = 1e-06, 
                            "delta_out"          = 1e-06,
                            "min_mu_k"           = 1e-06,
                            "min_Psi_d"          = 1e-06, 
                            "min_increase"       = 1e-04, 
                            #"norm_beta"          = FALSE, 
                            #"use_penalty"        = TRUE,
                            #"deriv_type"         = "theoretic",
                            #"step_halving"       = TRUE,
                            "delta_halving"      = 1e-06,
                            "min_iter_halving_in"  = 1,
                            "min_iter_halving_out" = 2,
                            "opt_function"       = "optim",
                            #"start_logit"        = TRUE,
                            #"start_grid_lambda"  = TRUE,
                            "initial_sm"         = FALSE,
                            "fix_smooth"         = FALSE, 
                            "sm_par_vec"         = c("lambda"=10,"s(x1)"=200,"s(x2)"=1000,"s(x3)"=2000,"s(x4)"=2000),
                            "sp_range"           = c(1e-8, 1e15),
                            #"max_iter_out_opt" = 25,
                            #"max_iter_in_opt"  = 25,
                            "reltol_opt"       = 1e-06,
                            #"less_warnings"      = TRUE,
                            "quietly"            = FALSE,
                            "save_step_response" = FALSE,
                            #"speed_up "          = FALSE,
                            "initial_model"      = "with_intercept")
        control_new$sm_par_vec           <- match.arg_sm(argu = control$sm_par_vec           , default_value = 10, formula = formula, data=data, initial_sm = control_new$initial_sm, fix_smooth=control_new$fix_smooth)
        
    } else {
        control_new <- list()
        if(!is.list(control)) { 
            stop("control must be a list") 
        } else {
            if(!all(names(control) %in% c("max_iter_in",
                                          "max_iter_out",
                                          "delta_in",
                                          "delta_out",
                                          "min_mu_k",
                                          "min_Psi_d",
                                          "min_increase",
                                          #"use_penalty",
                                          #"deriv_type",
                                          #"step_halving",
                                          "delta_halving",
                                          "min_iter_halving_in",
                                          "min_iter_halving_out",
                                          "opt_function",
                                          #"start_logit",
                                          #"start_grid_lambda",
                                          "initial_sm",
                                          "fix_smooth",
                                          "sm_par_vec",
                                          "sp_range",
                                          #"max_iter_out_opt",
                                          #"max_iter_in_opt",
                                          "reltol_opt",
                                         # "less_warnings",
                                          "quietly",
                                          "save_step_response",
                                          "initial_model"))) {
                stop("Unknown control parameter")
            }
            control_new$max_iter_in          <- match.arg_my(argu = control[["max_iter_in"]]          , default_value = 100)
            control_new$max_iter_out         <- match.arg_my(argu = control[["max_iter_out"]]         , default_value = 100)
            control_new$delta_in             <- match.arg_my(argu = control[["delta_in"]]             , default_value = 1e-06)
            control_new$delta_out            <- match.arg_my(argu = control[["delta_out"]]            , default_value = 1e-06)
            control_new$min_mu_k             <- match.arg_my(argu = control[["min_mu_k"]]             , default_value = 1e-06)
            control_new$min_Psi_d            <- match.arg_my(argu = control[["min_Psi_d"]]            , default_value = 1e-06)
            control_new$min_increase         <- match.arg_my(argu = control[["min_increase"]]         , default_value = 1e-04)
            #control_new$norm_beta            <- match.arg_my(argu = control[["norm_beta"]]            , default_value = FALSE)
            #control_new$use_penalty          <- match.arg_my(argu = control[["use_penalty"]]          , default_value = TRUE)
            #control_new$deriv_type           <- match.arg(control[["deriv_type"]] , choices = c("numeric","theoretic"))
            #control_new$step_halving         <- match.arg_my(argu = control[["step_halving"]]         , default_value = TRUE)
            control_new$delta_halving        <- match.arg_my(argu = control[["delta_halving"]]        , default_value = 1e-06)
            control_new$min_iter_halving_in  <- match.arg_my(argu = control[["min_iter_halving_in"]]  , default_value = 1)
            control_new$min_iter_halving_out <- match.arg_my(argu = control[["min_iter_halving_out"]] , default_value = 2)
            control_new$opt_function         <- match.arg(control[["opt_function"]] , choices = c("optim","nlminb"))
            control_new$initial_sm           <- match.arg_my(argu = control[["initial_sm"]]           , default_value =  TRUE)
            #control_new$start_logit          <- match.arg_my(argu = control[["start_logit"]]          , default_value =  TRUE)
            #control_new$start_grid_lambda    <- match.arg_my(argu = control[["start_grid_lambda"]]    , default_value =  TRUE)
            control_new$fix_smooth           <- match.arg_my(argu = control[["fix_smooth"]]           , default_value =  FALSE) 
            control_new$sm_par_vec           <- match.arg_sm(argu = control[["sm_par_vec"]]           , default_value = 10, formula = formula, data=data, initial_sm = control_new$initial_sm, fix_smooth=control_new$fix_smooth)
            control_new$sp_range             <- match.arg_ra(argu = control[["sp_range"]]             , default_value = c(1e-8, 1e15))
            #control_new$max_iter_out_opt   <- match.arg_my(argu = control[["max_iter_out_opt"]]   , default_value =  25)
            #control_new$max_iter_in_opt    <- match.arg_my(argu = control[["max_iter_in_opt"]]    , default_value =  25)
            control_new$reltol_opt         <- match.arg_my(argu = control[["reltol_opt"]]         , default_value =  1e-06)
            #control_new$less_warnings        <- match.arg_my(argu = control[["less_warnings"]]        , default_value =  TRUE)
            control_new$quietly              <- match.arg_my(argu = control[["quietly"]]              , default_value =  FALSE)
            control_new$save_step_response   <- match.arg_my(argu = control[["save_step_response"]]   , default_value =  FALSE)
            #control_new$speed_up             <- match.arg_my(argu = control[["speed_up"]]             , default_value =  TRUE)
            control_new$initial_model        <- match.arg(control[["initial_model"]] , choices = c("with_intercept","no_intercept"))
            
        }
        
        
    }
    control_new
}    

match.arg_my <- function(argu=NULL, default_value = NULL) {
    if(is.null(argu)){
        argu <- default_value
    } else {
        if(is.numeric(default_value)) {
            if(!is.numeric(argu)) {
                stop(paste( deparse(substitute(argu)), " must be numeric!"))
            } else {
                if(is.matrix(argu)) {
                    stop(paste( deparse(substitute(argu)), " must be a scalar!"))
                } else {
                    if(length(argu) > 1) {
                        stop(paste( deparse(substitute(argu)), " must be a scalar!"))
                    } else {
                        if(argu <= 0 ) { 
                            stop(deparse(substitute(argu)), " must be positive!")
                        } else {
                            argu <- argu
                        }
                    }
                }
            }
        } else {
            if(!is.logical(argu)) {
                stop(paste( deparse(substitute(argu)), " must be logical!"))
            } else {
                if(is.matrix(argu)) {
                    stop(paste( deparse(substitute(argu)), " must be a scalar!"))
                } else {
                    if(length(argu) > 1) {
                        stop(paste( deparse(substitute(argu)), " must be a scalar!"))
                    } else {
                        if(is.na(argu) ) {
                            stop(paste( deparse(substitute(argu)), " must be TRUE or FALSE!"))
                        }
                        argu <- argu
                    }
                }
            }
            
        }
    }
    argu 
}


match.arg_ra <- function(argu=NULL, default_value = NULL) {
    if(is.null(argu)) {
        argu <- default_value 
    } else {
        if(is.vector(argu, mode = "numeric")) {
            if(length(argu) != 2) {
                stop("range must be a vector of length 2")
            } else {
                if(argu[1] >= argu[2]) {
                    stop("range[1] must be smaller than range[2]")
                }
                if(any(argu <= 0)) {
                    stop("range must be positive")
                    }
                argu <- argu
            }
        } else {
            stop("range must be a numeric vector")
        }
    }
    argu
}

match.arg_sm <- function(argu = NULL, default_value = 10, formula = NULL, data=data,
                         initial_sm  = TRUE, fix_smooth = FALSE) {
    gp <- interpret.gam(formula) # interpret the formula
    gp_model <- gam(formula, data=data, fit=F)
    if(!is.null(argu)) {
        argu_short <- argu[-1]
        if(!("lambda" %in% names(argu))) {
            stop("lambda missing in sm_par_vec!")
        } else {
            if(!is.numeric(argu[which(names(argu) == "lambda")]) ) {
                stop(paste("lambda must be numeric!"))
            } else {
                if(argu[which(names(argu) == "lambda")] <= 0 ) {
                    stop(paste("lambda must be positive!"))
                }
            }
        }
        if(length(gp_model$sp) > 0){
            if(length(names(gp_model$sp)) != length(names(argu_short))) {
                stop("Wrong number of smoothing parameters")
            }
            if(!all(names(gp_model$sp) %in% names(argu_short))) {
                stop("Missing initial smoothing parameter!")
            } else {
                if(!all(names(gp_model$sp) == names(argu_short))) {
                    stop("Wrong ordering of smoothing parameters")
                }
            }
            if(!all(is.numeric(argu_short))) {
                stop("Initial smoothing parameter must be numeric!")
            }
            if(!all(argu_short > 0)) {
                stop("Initial smoothing parameter must be positive!")
            }
        } else {
            argu <- argu["lambda"]
        }
        
    } else {
        if(!initial_sm || fix_smooth) {
            warning(paste("Missing initial smoothing parameter values! They are set to ", default_value))
        }
        gp <- interpret.gam(formula) # interpret the formula 
        argu <- rep(default_value, times=length(gp_model$sp) + 1)
        names(argu) <- c("lambda",names(gp_model$sp))
    }
    argu   
}

