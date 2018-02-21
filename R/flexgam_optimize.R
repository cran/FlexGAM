flexgam_optimize_temp <- function(sm_par_log, formula_1, data_1, 
                                  #split_formula = split_formula, 
                                  split_list    = split_list,
                                  type          = type,
                                  family,
                                  try_error     = FALSE,
                                  lambda_default = 10,
                                  control_Opt = control_Opt)  {  
    
    
    control_Opt$sm_par_vec <- exp(sm_par_log)
    
    if(all(exp(sm_par_log) <= control_Opt$sp_range[2]) && all(exp(sm_par_log) >= control_Opt$sp_range[1])) {
        if(try_error) control_Opt$sm_par_vec["lambda"] <- lambda_default
        control_Opt$quietly    <- TRUE
        control_Opt$save_step_response <- FALSE
        
        #if(length(control_Opt$sm_par_vec) == 1) control_Opt$sm_par_vec <- c("lambda" = control_Opt$sm_par_vec)
        
        formula_2 <- formula_1
        UBRE_vec <- rep(NA,times=length(split_list))
        for(icv in 1:length(split_list)) {
            data_train <- data_1[-split_list[[icv]],]
            data_valid <- data_1[ split_list[[icv]],]
            
            G <- try(flexgam_inner(formula=formula_2, data=data_train, type=type, family = family, split_list = split_list, estimate_se = FALSE, 
                                   control = control_Opt),TRUE)
            if(!inherits(G, "try-error")) {
                UBRE_vec[icv] <- deviance.flexgam(object=G, newdata=data_valid)
                
            } else {
                UBRE_vec[icv] <- NA
            }
        }
        if(!all(is.na(UBRE_vec))){
            UBRE <- mean(UBRE_vec,na.rm=TRUE)
        } else {
            UBRE <- 1e10
            warning("All CV NA")
        }
    } else {
        UBRE <- 1e10
    }
    UBRE
}

flexgam_optimize <- function(formula_1a, data_1a, 
                             type        = c("FlexGAM1","FlexGAM2","FlexGAM1n","FlexGAM2n"),
                             family, 
                             split_list = split_list, 
                             control = control) {
    
    opt_function <- match.arg(control$opt_function,c("optim","nlminb"))
    
    control_in            <- control
    control_in$type       <- type
    control_in$formula    <- formula
    
    control_Opt <- control_in
    #control_Opt$less_warnings <- TRUE
    
    control_Opt$fix_smooth    <- TRUE
    
    if(is.null(split_list)) {
        split_list <- chunk3(vec=1:nrow(data_1a),ng = 5)
    }
    
    fixed <- FALSE
    trace_temp <- 0
    if(!control_in$quietly) {
        if(opt_function == "nlminb") trace_temp <- 3
        if(opt_function == "optim") trace_temp <- 1
    }
    try_error <- FALSE
    if(!fixed) {
        sm_par_log <- log(control_Opt$sm_par_vec)
        if(opt_function == "nlminb") {
            acv.min <- try(
                nlminb(objective = flexgam_optimize_temp, start = sm_par_log, 
                       formula_1 = formula_1a, 
                       data_1 = data_1a, 
                       split_list = split_list,
                       type = type,
                       family = family,
                       try_error = try_error,
                       lambda_default = control_Opt$sm_par_vec["lambda"],
                       control_Opt = control_Opt, 
                       lower = log(control_Opt$sp_range[1]), upper = log(control_Opt$sp_range[2]), 
                       control = list("trace" = trace_temp, rel.tol = control_Opt$reltol_opt))
                , TRUE)
            if(inherits(acv.min, "try-error")) {
                #print("B")
                try_error <- TRUE
                acv.min <- try(
                    nlminb(objective = flexgam_optimize_temp, start = sm_par_log[-"lambda"], 
                           formula_1 = formula_1a, 
                           data_1 = data_1a, 
                           split_list = split_list,
                           type = type,
                           family = family,
                           try_error = try_error,
                           lambda_default = control_Opt$sm_par_vec["lambda"],
                           control_Opt = control_Opt, 
                           lower = log(control_Opt$sp_range[1]), upper = log(control_Opt$sp_range[2]), 
                           control = list("trace" = trace_temp, rel.tol = control_Opt$reltol_opt))
                    , TRUE)
                if(inherits(acv.min, "try-error")) stop("Error in nlminb with fixed lambda!")
                
            }
            if(!try_error) {
                sp_vec_fin <- exp(acv.min$par)
            } else {
                sp_vec_fin <- c(control_Opt$sm_par_vec["lambda"],exp(acv.min$par))
            }
        }    
        if(opt_function == "optim") {
            #if(length(sm_par_log) > 1) {
            acv.min <- try(
                optim(fn = flexgam_optimize_temp, par = sm_par_log, formula_1 = formula_1a, 
                      data_1 = data_1a,
                      split_list = split_list,
                      type = type,
                      family = family,
                      try_error = try_error,
                      lambda_default = control_Opt$sm_par_vec["lambda"],
                      control_Opt = control_Opt, 
                      method="Nelder-Mead", control = list("reltol" = control_Opt$reltol_opt, "trace"=trace_temp))
                ,TRUE)
            if(inherits(acv.min, "try-error")) {
                try_error <- TRUE
                acv.min <- try(
                    optim(fn = flexgam_optimize_temp, par = sm_par_log[-"lambda"], formula_1 = formula_1a, 
                          data_1 = data_1a,
                          split_list = split_list,
                          type = type,
                          family = family,
                          try_error = try_error,
                          lambda_default = control_Opt$sm_par_vec["lambda"],
                          control_Opt = control_Opt, 
                          method="Nelder-Mead", control = list("reltol" = control_Opt$reltol_opt, "trace"=trace_temp))
                    ,TRUE)
                
                if(inherits(acv.min, "try-error")) stop("Error in optim with fixed lambda!")
            }
            if(!try_error) {
                sp_vec_fin <- exp(acv.min$par)
            } else {
                sp_vec_fin <- c(control_Opt$sm_par_vec["lambda"],exp(acv.min$par))
            }
        } #else {
        #     acv.min <- try(
        #         optim(sm_par_log, fn = flexgam_optimize_temp,  formula_1 = formula_1a, 
        #               data_1 = data_1a,
        #               split_list = split_list,
        #               type = type,
        #               family = family,
        #               try_error = try_error,
        #               lambda_default = control_Opt$sm_par_vec["lambda"],
        #               control_Opt = control_Opt, 
        #               lower = log(control_Opt$sp_range[1]), upper = log(control_Opt$sp_range[2]), 
        #               method="Brent", control = list("reltol" = control_Opt$reltol_opt, "trace"=trace_temp))
        #         ,TRUE)
        #     if(inherits(acv.min, "try-error")) {
        #         try_error <- TRUE
        #         acv.min <- try(
        #             optim(fn = flexgam_optimize_temp, par = sm_par_log[-"lambda"], formula_1 = formula_1a, 
        #                   data_1 = data_1a,
        #                   split_list = split_list,
        #                   type = type,
        #                   family = family,
        #                   try_error = try_error,
        #                   lambda_default = control_Opt$sm_par_vec["lambda"],
        #                   control_Opt = control_Opt, 
        #                   lower = log(control_Opt$sp_range[1]), upper = log(control_Opt$sp_range[2]), 
        #                   method="Brent", control = list("reltol" = control_Opt$reltol_opt, "trace"=trace_temp))
        #             ,TRUE)
        #         
        #         if(inherits(acv.min, "try-error")) stop("Error in optim with fixed lambda!")
        #     }
        #     if(!try_error) {
        #         sp_vec_fin <- c("lambda"=exp(acv.min$par))
        #     } else {
        #         sp_vec_fin <- c(control_Opt$sm_par_vec["lambda"])
        #     }
        #}
        #}
    } else {
        sp_vec_fin <- control_in$sm_par_vec
        acv.min <- NULL
    }
    
    formula_2 <- formula_1a
    control_Fin              <- control_in
    control_Fin$sm_par_vec   <- sp_vec_fin
    control_Fin$fix_smooth   <- TRUE
    
    if(!control_Fin$quietly){
        cat("\n","Final smoothing parameters jointly:\n")
        print(control_Fin$sm_par_vec)
    }    
    
    G_fin <- try(flexgam_inner(formula=formula_2, data=data_1a, type=type, family = family, split_list = split_list, estimate_se = TRUE,
                               control = control_Fin),TRUE)
    if(inherits(G_fin, "try-error")) {
        warning("G_fin with error refit with lambda=min(sp_range)")
        control_Fin$sm_par_vec["lambda"] <- control_Fin$sp_range[1]
        G_fin <- try(flexgam_inner(formula=formula_2, data=data_1a, type=type, family = family, split_list = split_list, estimate_se = TRUE,
                                   control = control_Fin),TRUE)
        if(inherits(G_fin, "try-error")) {
            warning("G_fin with error refit with lambda=min(sp_range)+0.005")
            control_Fin$sm_par_vec["lambda"] <- control_Fin$sp_range[1] + 0.005
            G_fin <- flexgam_inner(formula=formula_2, data=data_1a, type=type, family = family, split_list = split_list, estimate_se = TRUE,
                                   control = control_Fin)
        }
    }
    
    if(opt_function=="optim")   G_fin$CV_score = acv.min$value
    if(opt_function=="nlminb")  G_fin$CV_score = acv.min$objective
    G_fin
}


