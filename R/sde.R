
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
        #' @param response Vector of names of response variables
        #' @param fixpar Vector of names of fixed SDE parameters
        #' 
        #' @return A new SDE object
        initialize = function(formulas, data, type, response, 
                              par0 = NULL, fixpar = NULL) {
            private$formulas_ <- formulas
            private$type_ <- type
            private$response_ <- response
            private$fixpar_ <- fixpar
            
            if(any(!response %in% colnames(data)))
                stop("'response' not found in 'data'")
            
            # Link functions for SDE parameters
            link <- switch (type,
                            "BM" = list(mu = identity, sigma = log),
                            "OU" = list(mu = identity, beta = log, sigma = log),
                            "CTCRW" = list(beta = log, sigma = log),
                            "ESEAL_SSM" = list(mu = identity, sigma = log))
            private$link_ <- link
            
            # Inverse link functions for SDE parameters
            invlink <- switch (type,
                               "BM" = list(mu = identity, sigma = exp),
                               "OU" = list(mu = identity, beta = exp, sigma = exp),
                               "CTCRW" = list(beta = exp, sigma = exp),
                               "ESEAL_SSM" = list(mu = identity, sigma = exp))
            private$invlink_ <- invlink
            
            # Check that "formulas" is of the right length and with right names
            if(length(formulas) != length(invlink)) {
                err <- paste0("'formulas' should be a list of length ", 
                              length(invlink), " for the model ", type,
                              ", with components ", 
                              paste(names(invlink), collapse = ", "))
                stop(err)
            } else if(any(names(formulas) != names(invlink))) {
                err <- paste0("'formulas' should be a list with components ", 
                              paste(names(invlink), collapse = ", "))
                stop(err)
            }
            
            # Check that data has an "ID" column, and that it's a factor
            if(!any(colnames(data) == "ID")) {
                warning(paste("No ID column found in 'data',",
                              "assuming same ID for all observations"))
                data$ID <- factor(1)
            } else {
                data$ID <- factor(data$ID)
            }
            
            # Check that data has a "time" column
            if(!any(colnames(data) == "time")) {
                stop("'data' should have a time column")
            }
            private$data_ <- data
            
            # Save terms of model formulas
            mats <- self$make_mat()
            ncol_fe <- mats$ncol_fe
            ncol_re <- mats$ncol_re
            private$terms_ <- list(ncol_fe = ncol_fe,
                                   ncol_re = ncol_re,
                                   names_fe = colnames(mats$X_fe),
                                   names_re_all = colnames(mats$X_re),
                                   names_re = names(ncol_re))
            
            # Initial parameters (zero if par0 not provided)
            self$update_coeff_fe(rep(0, sum(ncol_fe)))
            self$update_coeff_re(rep(0, sum(ncol_re)))
            self$update_lambda(rep(1, length(ncol_re)))
            
            # Set initial fixed coefficients if provided (par0)
            if(!is.null(par0)) {
                # Number of SDE parameters
                n_par <- length(self$formulas())
                
                if(length(par0) != n_par) {
                    stop("'par0' should be of length ", n_par,
                         " with one entry for each SDE parameter (",
                         paste0(names(self$formulas()), collapse = ", "), ")")
                }
                
                # First column of each X_fe for each SDE parameter
                i0 <- c(1, cumsum(ncol_fe)[-n_par] + 1)
                
                # Apply link to get parameters on working scale
                private$coeff_fe_[i0] <- sapply(1:n_par, function(i) {
                    self$link()[[i]](par0[i])
                })
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
        
        #' @description Name(s) of response variable(s)
        response = function() {return(private$response_)},
        
        #' @description Name(s) of fixed parameter(s)
        fixpar = function() {return(private$fixpar_)},
        
        #' @description Link functions
        link = function() {return(private$link_)},
        
        #' @description Inverse link functions
        invlink = function() {return(private$invlink_)},
        
        #' @description Fixed effect parameters
        coeff_fe = function() {return(private$coeff_fe_)},
        
        #' @description Random effect parameters
        coeff_re = function() {return(private$coeff_re_)},
        
        #' @description Smoothness parameters
        lambda = function() {return(private$lambda_)},
        
        #' @description Variance components of smooth terms
        #' 
        #' @details This function transforms the smoothness parameter of
        #' each smooth term into a standard deviation, given by 
        #' SD = 1/sqrt(lambda). It is particularly helpful to get the
        #' standard deviations of independent normal random effects.
        vcomp = function() {return(1/sqrt(private$lambda_))},
        
        #' @description Terms of model formulas
        terms = function() {return(private$terms_)},
        
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
        
        #' @description Data frame of observations (subset response
        #' variables out of full data frame)
        obs = function() {
            self$data()[, self$response(), drop = FALSE]
        },
        
        #' @description Print equation for this model
        eqn = function() {
            switch (self$type(),
                    "BM" = "dZ(t) = mu dt + sigma dW(t)",
                    "OU" = "dZ(t) = beta (mu - Z(t)) dt + sigma dW(t)",
                    "CTCRW" = paste0("dV(t) = - beta V(t) dt + sigma dW(t)\n", 
                                     "dZ(t) = V(t) dt"),
                    "ESEAL_SSM" = paste0("dL(t) = mu dt + sigma dW(t)\n", 
                                         "Z(i) ~ N(a1 + a2 L(i)/R(i), tau^2/h(i))"))
            
        },
        
        ##############
        ## Mutators ##
        ##############
        #' @description Update fixed effect coefficients
        #' 
        #' @param new_coeff New coefficient vector
        update_coeff_fe = function(new_coeff) {
            private$coeff_fe_ <- new_coeff
            names(private$coeff_fe_) <- self$terms()$names_fe
        },
        
        #' @description Update random effect coefficients
        #' 
        #' @param new_coeff New coefficient vector
        update_coeff_re = function(new_coeff) {
            private$coeff_re_ <- new_coeff
            names(private$coeff_re_) <- self$terms()$names_re_all
        },
        
        #' @description Update smoothness parameters
        #' 
        #' @param new_coeff New smoothness parameter vector
        update_lambda = function(new_lambda) {
            private$lambda_ <- new_lambda
            names(private$lambda_) <- self$terms()$names_re
        },
        
        ###################
        ## Other methods ##
        ###################
        #' @description Print SDE and parameter formulas
        message = function() {
            message("#######################")
            message("### smoothSDE model ###")
            message("#######################")
            
            # Print SDE
            eqn <- self$eqn()
            message("SDE for ", self$type(), " model:")
            message(eqn, "\n")
            
            # Print parameter formulas
            message("Formulas for model parameters:")
            f <- self$formulas()
            for(i in 1:length(f)) {
                if(names(f)[i] %in% self$fixpar()) {
                    this_form <- "fixed"
                } else {
                    this_form <- as.character(f[[i]])[2]
                    this_form <- gsub("\\+", "+\n\t", this_form)
                }
                message(names(f)[i], " ~ ", this_form)
            }
            cat("\n")
        },
        
        #' @description Indices of fixed coefficients in coff_fe
        ind_fixcoeff = function() {
            # Number of SDE parameters
            n_par <- length(self$formulas())
            
            # Number of columns of X_fe for each SDE parameter
            ncol_fe <- self$make_mat()$ncol_fe
            
            # Counter for coefficients in coeff_fe
            k <- 1
            # Initialise vector of indices of fixed coefficients
            ind_fixcoeff <- NULL
            # Loop over SDE parameters
            for(par in 1:n_par) {
                # If this parameter is fixed, add corresponding indices
                # to ind_fixcoeff
                if(names(self$formulas())[par] %in% self$fixpar()) {
                    ind_thispar <- k:(k + ncol_fe[par] - 1)
                    ind_fixcoeff <- c(ind_fixcoeff, ind_thispar)
                }
                k <- k + ncol_fe[par]
            }
            
            return(ind_fixcoeff)
        },
        
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
        #'   \item ncol_fe Number of columns for X_fe for each parameter
        #'   \item ncol_re Number of columns of X_re and S for each random effect
        #' }
        make_mat = function(new_data = NULL) {
            # Initialise lists of matrices
            X_list_fe <- list()
            X_list_re <- list()
            S_list <- list()
            ncol_fe <- NULL
            ncol_re <- NULL
            names_fe <- NULL
            names_re <- NULL
            names_ncol_re <- NULL
            k <- 1
            
            # Loop over formulas
            for(j in seq_along(self$formulas())) {
                form <- self$formulas()[[j]]
                par_name <- names(self$formulas())[j]
                
                # Create matrices based on this formula
                if(is.null(new_data)) {
                    gam_setup <- gam(formula = update(form, dummy ~ .), 
                                     data = cbind(dummy = 1, self$data()), 
                                     fit = FALSE)
                    Xmat <- gam_setup$X
                    # Extract column names for design matrices
                    term_names <- gam_setup$term.names
                } else {
                    # Get design matrix for new data set
                    gam_setup <- gam(formula = update(form, dummy ~ .), 
                                     data = cbind(dummy = 1, self$data()))
                    Xmat <- predict(gam_setup, newdata = new_data, type = "lpmatrix")
                    # Extract column names for design matrices
                    term_names <- names(gam_setup$coefficients)
                }
                
                # Fixed effects design matrix
                X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
                subnames_fe <- paste0(par_name, ".", term_names[1:gam_setup$nsdf])
                names_fe <- c(names_fe, subnames_fe)
                
                # Random effects design matrix
                X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
                if(ncol(X_list_re[[k]]) > 0) {
                    subnames_re <- paste0(par_name, ".", term_names[-(1:gam_setup$nsdf)])
                    names_re <- c(names_re, subnames_re)                    
                }
                
                # Smoothing matrix
                S_list[[k]] <- bdiag_check(gam_setup$S)
                
                # Number of columns for fixed effects
                ncol_fe <- c(ncol_fe, gam_setup$nsdf)
                
                if(length(gam_setup$S) > 0) {
                    # Number of columns for each random effect
                    sub_ncol_re <- sapply(gam_setup$S, ncol)
                    ncol_re <- c(ncol_re, sub_ncol_re)
                    # Hacky way to get the names of smooth terms 
                    # (one for each column of ncol_re)
                    # regex from datascience.stackexchange.com/questions/8922
                    s_terms_i1 <- cumsum(c(1, sub_ncol_re[-length(sub_ncol_re)]))
                    s_terms <- gsub("(.*)\\..*", "\\1", subnames_re[s_terms_i1])
                    names_ncol_re <- c(names_ncol_re, s_terms)
                }
                
                k <- k + 1
            }
            
            # Store as block diagonal matrices
            X_fe <- bdiag_check(X_list_fe)
            colnames(X_fe) <- names_fe
            X_re <- bdiag_check(X_list_re)
            colnames(X_re) <- names_re
            S <- bdiag_check(S_list)
            
            # Name elements of ncol_re
            names(ncol_re) <- names_ncol_re
            
            return(list(X_fe = X_fe, X_re = X_re, S = S, 
                        ncol_fe = ncol_fe, ncol_re = ncol_re))
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
            # Number of time steps
            n <- nrow(self$data())
            
            # Create model matrices
            mats <- self$make_mat()
            X_fe <- mats$X_fe
            X_re <- mats$X_re
            S <- mats$S
            ncol_fe <- mats$ncol_fe
            ncol_re <- mats$ncol_re
            
            # Format initial parameters for TMB
            # (First fixed effects, then random effects)
            tmb_par <- list(coeff_fe = self$coeff_fe(),
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
                tmb_par$coeff_re <- self$coeff_re()
                tmb_par$log_lambda <- rep(0, length(ncol_re))
            }
            
            # Setup fixed parameters
            if(!is.null(self$fixpar())) {
                # Indices of fixed coefficients in coeff_fe
                ind_fixcoeff <- self$ind_fixcoeff()
                
                # Define vector with a different integer for each coefficient
                # to be estimated, and NA for each fixed coefficient
                coeff_fe_map <- 1:ncol(X_fe)
                coeff_fe_map[ind_fixcoeff] <- NA
                
                # Update map (to be passed to TMB)
                map <- c(map, list(coeff_fe = factor(coeff_fe_map)))
            }
            
            # TMB data object
            tmb_dat <- list(type = self$type(),
                            ID = self$data()$ID,
                            times = self$data()$time,
                            obs = as.matrix(self$obs()),
                            X_fe = X_fe,
                            X_re = X_re,
                            S = S,
                            ncol_re = ncol_re)
            
            # Define initial state and covariance for Kalman filter
            if(self$type() == "CTCRW") {
                # First index for each ID
                i0 <- c(1, which(self$data()$ID[-n] != self$data()$ID[-1]) + 1)
                # Initial state = (x1, 0, y1, 0)
                tmb_dat$a0 <- cbind(self$obs()[i0, 1], 0, self$obs()[i0,2], 0)
                # Initial state covariance
                tmb_dat$P0 <- diag(c(1, 10, 1, 10))
            } else if(self$type() == "ESEAL_SSM") {
                # Initial state = initial lipid mass, obtained from
                # Depart_Arrive_FatMass_measurements.csv in the supporting
                # information of Pirotta et al. (2019, Behav. Ecol.)
                dep_fat <- c(34, 40.3, 38.6, 31.5, 40.4, 40.2, 29.9, 41.8, 28.3, 40.8,
                             44, 35.2, 44.9, 40.2, 32.1, 38.3, 35.2, 39.5, 35.6, 30.4,
                             41, 32.9, 41.6, 47.1, 34.8, 36.1)
                tmb_dat$a0 <- cbind(1, dep_fat)
                tmb_dat$P0 <- diag(c(0, 10))
                
                # Initialise model-specific parameters
                ssm_par <- list(log_tau = log(1),
                                a1 = -0.578,
                                log_a2 = log(1.214))
                tmb_par <- c(ssm_par, tmb_par)
                
                # Number of daily drift dives
                tmb_dat$h <- self$data()$h
                # Non-lipid tissue mass
                tmb_dat$R <- self$data()$R
            }
            
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
            # Print model formulation
            self$message()
            
            # Setup if necessary
            if(is.null(private$tmb_obj_)) {
                self$setup(silent = silent)
            }
            
            # Fit model
            private$fit_ <- do.call(optim, private$tmb_obj_)
            # Get estimates and precision matrix for all parameters
            private$tmb_rep_ <- sdreport(private$tmb_obj_, getJointPrecision = TRUE, 
                                         skip.delta.method = FALSE)
            
            # Save parameter estimates
            par_list <- as.list(private$tmb_rep_, "Estimate")
            self$update_coeff_fe(par_list$coeff_fe)
            if(length(self$terms()$ncol_re) > 0) {
                # Only save coeff_re and lambda is there are random effects
                self$update_coeff_re(par_list$coeff_re)
                self$update_lambda(exp(par_list$log_lambda))
            }
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
            # Number of SDE parameters
            n_par <- length(self$formulas())
            # Number of time steps
            n <- nrow(X_fe)/n_par
            # Indices of non-fixed SDE parameters
            ind_estpar <- which(!names(self$formulas()) %in% self$fixpar())
            # Indices of estimated coefficients in coeff_fe
            fe_cols <- rep(1:n_par, sde$terms()$ncol_fe)
            ind_est_fe <- which(fe_cols %in% ind_estpar)
            
            # TMB report
            rep <- self$tmb_rep()
            
            # Joint covariance matrix
            if(!is.null(rep$jointPrecision)) {
                jointCov <- as.matrix(solve(rep$jointPrecision))
                colnames(jointCov) <- colnames(rep$jointPrecision)
                rownames(jointCov) <- colnames(jointCov)
            } else {
                # If there are no random effects
                jointCov <- rep$cov.fixed
            }
            
            # Vector of all parameters
            par_all <- c(rep$par.fixed, rep$par.random)
            
            # Make sure that parameters have the same order in vector of estimates
            # and in covariance matrix
            if(!all(names(par_all) == colnames(jointCov)))
                stop("Check TMB parameter order (should be fixed first, then random)")
            
            # Posterior draws from MVN(par_all, jointCov)
            par_post <- rmvn(n = n_post, mu = par_all, V = jointCov)
            
            # In post_fe, set columns for fixed parameters to fixed value,
            # and use posterior draws for non-fixed parameters
            post_fe <- matrix(rep(self$coeff_fe(), each = n_post), 
                              nrow = n_post, ncol = sum(sde$terms()$ncol_fe))
            post_fe[,ind_est_fe] <- par_post[, which(colnames(par_post) == "coeff_fe")]
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
        
        #' @description Model residuals
        residuals = function() {
            data <- self$data()
            n <- nrow(data)
            
            # Get start and end indices of tracks
            break_ind <- which(data$ID[-1] != data$ID[-n])
            start_ind <- c(1, break_ind + 1)
            end_ind <- c(break_ind, n)
            
            # Time intervals
            dtimes <- data$time[-start_ind] - data$time[-end_ind]
            
            # Get SDE parameters for each time step
            mats <- self$make_mat()
            par <- self$par_all(X_fe = mats$X_fe, X_re = mats$X_re)
            
            # Response variable
            Z <- data[[self$response()]]
            
            # Compute mean and sd of normal transition density
            if(self$type() == "BM") {
                mean <- Z[-end_ind] + par[-end_ind, "mu"] * dtimes
                sd <- par[-end_ind, "sigma"] * sqrt(dtimes)
            } else if(self$type() == "OU") {
                mu <- par[-end_ind, "mu"]
                beta <- par[-end_ind, "beta"]
                sigma <- par[-end_ind, "sigma"]
                mean <-  mu + exp(- beta * dtimes) * (Z[-end_ind] - mu)
                sd <- sigma/sqrt(2 * beta) * sqrt(1 - exp(-2 * beta * dtimes))
            } else {
                stop(paste("Residuals not implemented for model", self$type()))
            }
            
            # Residuals ~ N(0, 1) under assumptions of model and discretization
            res <- (Z[-start_ind] - mean) / sd
            return(res)
        },
        
        #' @description Akaike Information Criterion
        #'
        #' This function is adapted from Dave Miller's code in the
        #' package CTMCdive
        #'
        #' @param k Penalty per parameter; the default 
        #' \code{k = 2} is the classical AIC
        AIC = function(k = 2) {
            llk <- - self$res()$value
            mats <- self$make_mat()
            
            # Get effective degrees of freedom for smooths
            lambda <- self$lambda()
            edf_re <- edf(X_re = as.matrix(mats$X_re), 
                          S = as.matrix(mats$S), 
                          lambda = lambda,
                          ncol_re = mats$ncol_re)
            
            # Degrees of freedom for fixed effects
            edf_fe <- length(self$res()$par) - length(lambda)
            
            # Total EDF
            edf_total <- edf_re + edf_fe
            
            AIC <- - 2 * llk + k * edf_total
            return(AIC)
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
                post <- self$post(X_fe = mats$X_fe, 
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
        response_ = NULL,
        fixpar_ = NULL,
        link_ = NULL,
        invlink_ = NULL,
        coeff_fe_ = NULL,
        coeff_re_ = NULL,
        lambda_ = NULL,
        terms_ = NULL,
        tmb_obj_ = NULL,
        fit_ = NULL,
        tmb_rep_ = NULL
    )
)
