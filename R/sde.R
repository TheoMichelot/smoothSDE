
#' R6 class for stochastic differential equation
#' 
#' Contains the model formulas and data.
SDE <- R6Class(
    classname = "SDE",
    
    public = list(
        #################
        ## Constructor ##
        #################
        #' @description Create a SDE object
        #' 
        #' @param formulas List of formulas for model parameters
        #' @param data Data frame with covariates, response variable,
        #' time, and ID
        #' @param type Type of SDE ("BM", "OU"...)
        #' 
        #' @return A new SDE object
        initialize = function(formulas, data, type) {
            private$formulas_ <- formulas
            private$data_ <- data
            private$type_ <- type
            
            # SDE type code (to pass to C++)
            type_code <- switch (type,
                                 "BM" = 1,
                                 "OU" = 2,
                                 stop("Invalid 'type'"))
            private$type_code_ <- type_code
            
            # Inverse link functions for SDE parameters
            invlink <- switch (type,
                               "BM" = list(mu = identity, sigma = exp),
                               "OU" = list(mu = identity, beta = exp, sigma = exp))
            private$invlink_ <- invlink
            
            # Check that "formulas" is of the right length
            if(length(formulas) != length(invlink)) {
                err <- paste0("'formulas' should be a list of length ", 
                              length(invlink), " for the model ", type,
                              ", with components ", 
                              paste(names(invlink), collapse = ", "))
                stop(err)
            }
        },
        
        ###############
        ## Accessors ##
        ###############
        #' @description Formulas of SDE object
        formulas = function() {return(private$formulas_)},
        
        #' @description Data of SDE object
        data = function() {return(private$data_)},
        
        #' @description Type of SDE object
        type = function() {return(private$type_)},
        
        #' @description Type of SDE object as integer code
        type_code = function() {return(private$type_code_)},
        
        #' @description Inverse link functions
        invlink = function() {return(private$invlink_)},
        
        #' @description Fixed effect parameters
        coeff_fe = function() {return(private$coeff_fe_)},
        
        #' @description Random effect parameters
        coeff_re = function() {return(private$coeff_re_)},
        
        #' @description Output of optimiser after model fitting
        res = function() {
            if (is.null(private$fit_)) {
                stop("Fit model first")
            }
            
            return(private$fit_)
        },
        
        #' @description Model object created by TMB. This is the output of 
        #' the TMB function \code{MakeADFun}, and it is a list including elements
        #' \itemize{
        #'   \item{\code{fn}}{Objective function}
        #'   \item{\code{gr}}{Gradient function of fn}
        #'   \item{\code{par}}{Vector of initial parameters on working scale}
        #' }
        tmb_obj = function() {
            if(is.null(private$tmb_obj_)) {
                stop("Setup model first")
            }
            
            return(private$tmb_obj_)
        },
        
        #' @description Output of the TMB function \code{sdreport}, which includes 
        #' estimates and standard errors for all model parameters.
        tmb_rep = function() {
            if(is.null(private$tmb_rep_)) {
                stop("Fit model first")
            }
            
            return(private$tmb_rep_)
        },
        
        ###################
        ## Other methods ##
        ###################
        #' @description Create model matrices
        #'
        #' @param new_data Optional new data set, including covariates for which
        #' the design matrices should be created.
        #' 
        #' @return A list of
        #' \itemize{
        #'   \item X_fe Design matrix for fixed effects
        #'   \item X_re Design matrix for random effects
        #'   \item S Smoothness matrix
        #'   \item ncol_re Number of columns of X_re and S for each random effect
        #' }
        make_mat = function(new_data = NULL) {
            # Initialise lists of matrices
            X_list_fe <- list()
            X_list_re <- list()
            S_list <- list()
            ncol_re <- NULL
            k <- 1
            
            # Loop over formulas
            for(form in self$formulas()) {
                # Create matrices based on this formula
                if(is.null(new_data)) {
                    gam_setup <- gam(formula = update(form, dummy ~ .), 
                                     data = cbind(dummy = 1, self$data()), 
                                     fit = FALSE)
                    Xmat <- gam_setup$X
                } else {
                    # Get design matrix for new data set
                    gam_setup <- gam(formula = update(form, dummy ~ .), 
                                     data = cbind(dummy = 1, self$data()))
                    Xmat <- predict(gam_setup, newdata = new_data, type = "lpmatrix")
                }
                
                # Fixed effects design matrix
                X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
                
                # Random effects design matrix
                X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
                
                # Smoothing matrix
                S_list[[k]] <- bdiag_check(gam_setup$S)
                
                # Number of columns for each random effect
                if(length(gam_setup$S) > 0)
                    ncol_re <- c(ncol_re, sapply(gam_setup$S, ncol))
                
                k <- k + 1
            }
            
            # Store as block diagonal matrices
            X_fe <- bdiag_check(X_list_fe)
            X_re <- bdiag_check(X_list_re)
            S <- bdiag_check(S_list)
            
            return(list(X_fe = X_fe, X_re = X_re, S = S, ncol_re = ncol_re))
        },
        
        #' Design matrices for grid of covariates
        #' 
        #' @param var Name of variable
        #' @param covs Optional data frame with a single row and one column
        #' for each covariate, giving the values that should be used. If this is
        #' not specified, the mean value is used for numeric variables, and the
        #' first level for factor variables.
        #' 
        #' @return A list with the same elements as the output of make_mat, 
        #' plus a data frame of covariates values.
        make_mat_grid = function(var, covs = NULL) {
            # Data frame for covariate grid
            new_data <- cov_grid(var = var, data = self$data(), covs = covs, 
                                 formulas = self$formulas())
            
            # Create design matrices
            mats <- self$make_mat(new_data = new_data)
            
            # Save data frame of covariate values
            mats$new_data <- new_data
            
            return(mats)
        },
        
        #' @description TMB setup
        #'  
        #' This creates an attribute \code{tmb_obj}, which can be used to 
        #' evaluate the negative log-likelihood function.
        #' 
        #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
        setup = function(silent = TRUE) {
            # Create model matrices
            mats <- self$make_mat()
            X_fe <- mats$X_fe
            X_re <- mats$X_re
            S <- mats$S
            ncol_re <- mats$ncol_re
            
            # Format initial parameters for TMB
            # (First fixed effects, then random effects)
            tmb_par <- list(coeff_fe = rep(0, ncol(X_fe)),
                            log_lambda = 0,
                            coeff_re = 0)
            
            # Setup random effects
            map <- NULL
            random <- NULL
            if(is.null(S)) {
                # If there are no random effects, 
                # coeff_re and log_lambda are not estimated
                map <- c(map, list(coeff_re = factor(NA),
                                   log_lambda = factor(NA)))
                S <- as(matrix(0, 1, 1), "sparseMatrix")
                ncol_re <- 0
                X_re <- as(rep(0, nrow(X_fe)), "sparseMatrix")
            } else {
                # If there are random effects, 
                # set initial values for coeff_re and log_lambda
                random <- c(random, "coeff_re")
                tmb_par$coeff_re <- rep(0, ncol(S))
                tmb_par$log_lambda <- rep(0, length(ncol_re))
            }
            
            # TMB data object
            tmb_dat <- list(ID = self$data()$ID,
                            times = self$data()$time,
                            Z = as.matrix(self$data()$Z),
                            X_fe = X_fe,
                            X_re = X_re,
                            S = S,
                            ncol_re = ncol_re,
                            type = self$type_code())
            
            # Create TMB object
            tmb_obj <- MakeADFun(data = tmb_dat, parameters = tmb_par, 
                                 dll = "smoothSDE", silent = silent,
                                 map = map, random = random)
            
            # Negative log-likelihood function
            private$tmb_obj_ <- tmb_obj
        },
        
        #' @description Model fitting
        #' 
        #' The negative log-likelihood of the model is minimised using the
        #' function \code{optim}. TMB uses the Laplace approximation to integrate 
        #' the random effects out of the likelihood.
        #' 
        #' After the model has been fitted, the output of \code{optim} can be
        #' accessed using the method \code{res}.
        #' 
        #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
        fit = function(silent = TRUE) {
            # Setup if necessary
            if(is.null(private$tmb_obj_)) {
                self$setup(silent = silent)
            }
            
            # Fit model
            private$fit_ <- do.call(optim, private$tmb_obj_)
            # Get estimates and precision matrix for all parameters
            private$tmb_rep_ <- sdreport(private$tmb_obj_, getJointPrecision = TRUE, 
                                         skip.delta.method = FALSE)
            
            # Save parameters
            par_list <- as.list(private$tmb_rep_, "Estimate")
            private$coeff_fe_ <- par_list$coeff_fe
            private$coeff_re_ <- par_list$coeff_re
        },
        
        #' @description Get parameters from design matrices
        #' 
        #' @param X_fe Design matrix for fixed effects, as returned
        #' by \code{make_mat}
        #' @param X_re Design matrix for random effects, as returned
        #' by \code{make_mat}
        #' 
        #' @return Matrix with one column for each parameter
        par_all = function(X_fe, X_re) {
            # Get linear predictor and put into matrix where each row
            # corresponds to a time step and each column to a parameter
            lp <- X_fe %*% self$coeff_fe() + X_re %*% self$coeff_re()
            lp_mat <- matrix(lp, ncol = length(self$formulas()))
            
            # Apply inverse link to get parameters on natural scale
            par_mat <- matrix(NA, nrow = nrow(lp_mat), ncol = ncol(lp_mat))
            for(i in 1:ncol(lp_mat)) {
                par_mat[,i] <- self$invlink()[[i]](lp_mat[,i])
            }
            colnames(par_mat) <- names(self$invlink())
            
            return(par_mat)
        },
        
        #' @description Posterior draws for uncertainty quantification
        #' 
        #' @param X_fe Design matrix (fixed effects)
        #' @param X_re Design matrix (random effects)
        #' @param n_post Number of posterior draws 
        post = function(X_fe, X_re, n_post = 100) {
            # Number of parameters
            n_par <- length(self$formulas())
            # Number of time steps
            n <- nrow(X_fe)/n_par
            
            # TMB report
            rep <- self$tmb_rep()
            
            # Joint covariance matrix
            jointCov <- as.matrix(solve(rep$jointPrecision))
            colnames(jointCov) <- colnames(rep$jointPrecision)
            
            # Vector of all parameters
            par_all <- c(rep$par.fixed, rep$par.random)
            
            # Make sure that parameters have the same order in vector of estimates
            # and in covariance matrix
            if(!all(names(par_all) == colnames(jointCov)))
                stop("Check TMB parameter order (should be fixed first, then random)")
            
            # Posterior draws from MVN(par_all, jointCov)
            par_post <- rmvn(n = n_post, mu = par_all, V = jointCov)
            
            post_fe <- par_post[, which(colnames(par_post) == "coeff_fe")]
            post_re <- par_post[, which(colnames(par_post) == "coeff_re")]
            lp_post <- X_fe %*% t(post_fe) + X_re %*% t(post_re)
            
            lp_array <- array(lp_post, dim = c(n, n_par, n_post))
            
            par_array <- array(NA, dim = c(n, n_par, n_post))
            for(i in 1:dim(lp_array)[2]) {
                par_array[,i,] <- self$invlink()[[i]](lp_array[,i,])
            }
            dimnames(par_array)[[2]] <- names(self$invlink())
            
            return(par_array)
        },
        
        ######################
        ## Plotting methods ##
        ######################
        #' @description Plot observation parameters
        #' 
        #' @param var Name of covariate as a function of which the parameters
        #' should be plotted
        #' @param covs Optional data frame with a single row and one column
        #' for each covariate, giving the values that should be used. If this is
        #' not specified, the mean value is used for numeric variables, and the
        #' first level for factor variables.
        #' @param n_post Number of posterior draws to plot. Default: 0.
        #' 
        #' @return A ggplot object
        plot_par = function(var, covs = NULL, n_post = 0) {
            # Create design matrices
            mats <- self$make_mat_grid(var = var, covs = covs)
            par <- self$par_all(X_fe = mats$X_fe, X_re = mats$X_re)
            
            # Data frame for posterior draws
            if(n_post > 0) {
                post <- my_sde$post(X_fe = mats$X_fe, 
                                    X_re = mats$X_re, 
                                    n_post = n_post)
                post_df <- as.data.frame.table(post)
                colnames(post_df) <- c("var", "par", "stratum", "val")
                
                post_df$mle = "no"
            } else {
                post_df <- NULL
            }
            
            # Data frame for MLE
            mle_df <- as.data.frame.table(par)
            colnames(mle_df) <- c("var", "par", "val")
            mle_df$stratum <- "mle"
            mle_df$mle <- "yes"
            
            # Full data frame
            df <- rbind(mle_df, post_df)
            df$var <- mats$new_data[, var]
            
            # Create caption with values of other (fixed) covariates      
            plot_txt <- NULL
            if(ncol(mats$new_data) > 1) {
                other_covs <- mats$new_data[1, which(colnames(mats$new_data) != var), 
                                            drop = FALSE]
                
                # Round numeric values, and transform factors to strings
                num_ind <- sapply(other_covs, is.numeric)
                other_covs[num_ind] <- lapply(other_covs[num_ind], function(cov) 
                    round(cov, 2))
                fac_ind <- sapply(other_covs, is.factor)
                other_covs[fac_ind] <- lapply(other_covs[fac_ind], as.character)
                
                plot_txt <- paste(colnames(other_covs), "=", other_covs, 
                                  collapse = ", ")
            }
            
            # Create plot
            pal <- c("no" = rgb(0.7, 0, 0, 0.1), "yes" = rgb(0, 0, 0, 1))
            p <- ggplot(df, aes(var, val, group = stratum, col = mle)) + 
                theme_light() + geom_line() +
                scale_colour_manual(values = pal, guide = "none") +
                facet_wrap(c("par"), scales = "free_y",
                           strip.position = "left",
                           labeller = label_bquote(.(as.character(par)))) +
                xlab(var) + ylab(NULL) + ggtitle(plot_txt) +
                theme(strip.background = element_blank(),
                      strip.placement = "outside", 
                      strip.text = element_text(colour = "black"))
            
            return(p)
        }
    ),
    
    private = list(
        formulas_ = NULL,
        data_ = NULL,
        type_ = NULL,
        type_code_ = NULL,
        invlink_ = NULL,
        coeff_fe_ = NULL,
        coeff_re_ = NULL,
        tmb_obj_ = NULL,
        fit_ = NULL,
        tmb_rep_ = NULL
    )
)
