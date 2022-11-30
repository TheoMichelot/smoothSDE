
#' R6 class for stochastic differential equation
#' 
#' Contains the model formulas and data.
#' 
#' @importFrom R6 R6Class
#' @importFrom mgcv gam rmvn
#' @importFrom ggplot2 ggplot aes theme_light geom_line theme scale_colour_manual
#' facet_wrap label_bquote xlab ylab ggtitle element_blank element_text geom_point
#' geom_ribbon scale_size_manual geom_histogram geom_vline
#' @importFrom TMB MakeADFun sdreport
#' 
#' @useDynLib smoothSDE, .registration = TRUE
#' 
#' @export
SDE <- R6Class(
    classname = "SDE",
    
    public = list(
        #################
        ## Constructor ##
        #################
        #' @description Create a SDE object
        #' 
        #' @param formulas List of formulas for model parameters, with one element
        #' for each SDE parameter. Formulas can use standard R syntax, as well
        #' as mgcv-style syntax for splines and random effects.
        #' @param data Data frame with covariates, response variable,
        #' time, and ID
        #' @param type Type of SDE. Options are "BM" (Brownian motion), "OU" (Ornstein-
        #' Uhlenbeck process), "CTCRW" (continuous-time correlated random walk, a.k.a.
        #' integrated Ornstein-Uhlenbeck process), "CIR" (Cox-Ingersoll-Ross process),
        #' "BM_SSM" (BM with measurement error), "OU_SSM" (OU with measurement error),
        #' "BM_t" (BM with Student's t-distributed increments)
        #' @param response Name of response variable, correspond to a column name in
        #' \code{data}. Can be a vector of names if multiple response variables
        #' @param par0 Vector of initial values for SDE parameters, with one value
        #' for each SDE parameter. If not provided, parameters are initialised to
        #' zero on the link scale.
        #' @param fixpar Vector of names of fixed SDE parameters
        #' @param other_data Named list of data objects to pass to likelihood, only
        #' required for special models
        #' 
        #' @return A new SDE object
        initialize = function(formulas = NULL, data, type, response, par0 = NULL, 
                              fixpar = NULL, other_data = NULL) {
            private$type_ <- type
            private$response_ <- response
            private$fixpar_ <- fixpar
            
            if(any(!response %in% colnames(data)))
                stop("'response' not found in 'data'")
            
            # Link functions for SDE parameters
            n_dim <- length(response)
            link <- switch (type,
                            "BM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                             sigma = log)),
                            "BM_SSM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                 sigma = log)),
                            "BM_t" = list(mu = identity, sigma = log),
                            "OU" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                             tau = log, kappa = log)),
                            "OU_SSM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                 tau = log, kappa = log)),
                            "CIR" = as.list(c(mu = lapply(1:n_dim, function(i) log), 
                                              beta = log, sigma = log)),
                            "CTCRW" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                tau = log, nu = log)),
                            "ESEAL_SSM" = list(mu = identity, sigma = log))
            
            # Inverse link functions for SDE parameters
            invlink <- switch (type,
                               "BM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                sigma = exp)),
                               "BM_SSM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                    sigma = exp)),
                               "BM_t" = list(mu = identity, sigma = exp),
                               "OU" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                tau = exp, kappa = exp)),
                               "OU_SSM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                    tau = exp, kappa = exp)),
                               "CIR" = as.list(c(mu = lapply(1:n_dim, function(i) exp), 
                                                 beta = exp, sigma = exp)),
                               "CTCRW" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                   tau = exp, nu = exp)),
                               "ESEAL_SSM" = list(mu = identity, sigma = exp))
            
            private$link_ <- link
            private$invlink_ <- invlink
            
            # Check that "formulas" is of the right length and with right names
            if(is.null(formulas)) {
                formulas <- lapply(invlink, function(f) return(~1))
            } else if(length(formulas) != length(invlink)) {
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
            if(any(sapply(formulas[fixpar], function(f) {f != ~1}))) {
                stop("formulas should be ~1 for fixed parameters")
            }
            private$formulas_ <- formulas
            
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
            
            # Save terms of model formulas and model matrices
            mats <- self$make_mat()
            ncol_fe <- mats$ncol_fe
            ncol_re <- mats$ncol_re
            private$terms_ <- list(ncol_fe = ncol_fe,
                                   ncol_re = ncol_re,
                                   names_fe = colnames(mats$X_fe),
                                   names_re_all = colnames(mats$X_re),
                                   names_re = names(ncol_re))
            private$mats_ <- list(X_fe = mats$X_fe, X_re = mats$X_re, S = mats$S)
            
            # Initial parameters (zero if par0 not provided)
            self$update_coeff_fe(rep(0, sum(ncol_fe)))
            self$update_coeff_re(rep(0, sum(ncol_re)))
            self$update_lambda(rep(1, length(ncol_re)))
            
            # Set initial fixed effect coefficients if provided (par0)
            if(!is.null(par0)) {
                # Number of SDE parameters
                n_par <- length(self$formulas())
                
                if(length(par0) != n_par) {
                    stop("'par0' should be of length ", n_par,
                         " with one entry for each SDE parameter (",
                         paste0(names(self$formulas()), collapse = ", "), ")")
                }
                
                # First column of X_fe for each SDE parameter
                i0 <- c(1, cumsum(ncol_fe)[-n_par] + 1)
                
                # Apply link to get parameters on working scale
                private$coeff_fe_[i0] <- sapply(1:n_par, function(i) {
                    self$link()[[i]](par0[i])
                })
            }
            
            # Process decay terms
            if(!is.null(other_data$t_decay)) {
                # Find columns for decay if necessary
                if(is.null(other_data$col_decay)) {
                    decay_term <- other_data$decay_term
                    str <- substr(self$terms()$names_re_all, 1, nchar(decay_term))
                    other_data$col_decay <- which(str == decay_term)                    
                }
                if(length(other_data$t_decay) != length(formulas) * nrow(data)) {
                    stop(paste0("'other_data$t_decay' should be of length (number ",
                                "of parameters) x (number of data)"))
                }
                if(length(other_data$col_decay) != length(other_data$ind_decay)) {
                    stop("Check length of 'other_data$ind_decay' and 'other_data$col_decay'")
                }
                private$rho_ <- rep(1, length(unique(other_data$ind_decay)))
            } else {
                private$rho_ <- 1
            }
            private$other_data_ <- other_data
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
        
        #' @description List of model matrices (X_fe, X_re, and S)
        mats = function() {return(private$mats_)},
        
        #' @description Named list of additional data objects
        other_data = function() {return(private$other_data_)},
        
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
        
        #' @description Standard deviations of smooth terms
        #' 
        #' This function transforms the smoothness parameter of
        #' each smooth term into a standard deviation, given by 
        #' SD = 1/sqrt(lambda). It is particularly helpful to get the
        #' standard deviations of independent normal random effects.
        sdev = function() {return(1/sqrt(private$lambda_))},
        
        #' @description Decay parameter
        rho = function() {
            return(private$rho_)
        },
        
        #' @description Terms of model formulas
        terms = function() {return(private$terms_)},
        
        #' @description Output of optimiser after model fitting
        out = function() {
            if (is.null(private$out_)) {
                stop("Fit model first")
            }
            
            return(private$out_)
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
        
        #' @description Model object created by TMB for the joint likelihood of
        #' the fixed and random effects. This is the output of the TMB function 
        #' \code{MakeADFun}, and it is a list including elements
        #' \itemize{
        #'   \item{\code{fn}}{Objective function}
        #'   \item{\code{gr}}{Gradient function of fn}
        #'   \item{\code{par}}{Vector of initial parameters on working scale}
        #' }
        tmb_obj_joint = function() {
            if(is.null(private$tmb_obj_joint_)) {
                stop("Setup model first")
            }
            
            return(private$tmb_obj_joint_)
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
        
        #' @description Get design matrix for random effects in decay model
        #' 
        #' The design matrix is obtained by taking X_re (returned by 
        #' make_mat), and multiplying the relevant columns by something like
        #' exp(-rho * time) to force the splines to decay to zero with a rate
        #' determined by rho.
        #' 
        #' @return Design matrix
        X_re_decay = function() {
            # Check that there are decaying terms
            if(!is.null(self$other_data()$t_decay)) {
                # Make design matrices
                X_re <- self$mats()$X_re
                
                # Decay parameters
                rho <- self$rho()
                t_decay <- self$other_data()$t_decay
                col_decay <- self$other_data()$col_decay
                ind_decay <- self$other_data()$ind_decay
                
                # Apply decay
                for(i in seq_along(col_decay)) {
                    col <- col_decay[i]
                    ind <- ind_decay[i]
                    X_re[,col] <- X_re[,col] * exp(-rho[ind] * t_decay)
                }
            } else {
                stop("This model has no decaying terms")
            }
            
            return(X_re)
        },
        
        ##############
        ## Mutators ##
        ##############
        #' @description Update fixed effect coefficients
        #' 
        #' @param new_coeff New coefficient vector
        update_coeff_fe = function(new_coeff) {
            private$coeff_fe_ <- matrix(new_coeff)
            rownames(private$coeff_fe_) <- self$terms()$names_fe
        },
        
        #' @description Update random effect coefficients
        #' 
        #' @param new_coeff New coefficient vector
        update_coeff_re = function(new_coeff) {
            private$coeff_re_ <- matrix(new_coeff)
            rownames(private$coeff_re_) <- self$terms()$names_re_all
        },
        
        #' @description Update smoothness parameters
        #' 
        #' @param new_lambda New smoothness parameter vector
        update_lambda = function(new_lambda) {
            private$lambda_ <- matrix(new_lambda)
            rownames(private$lambda_) <- self$terms()$names_re
        },
        
        #' @description Update decay parameter
        #' 
        #' @param new_rho New decay parameter vector
        update_rho = function(new_rho) {
            private$rho_ <- new_rho
        },
        
        #########################
        ## Make model matrices ##
        #########################
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
                X_list_fe[[j]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
                subnames_fe <- paste0(par_name, ".", term_names[1:gam_setup$nsdf])
                names_fe <- c(names_fe, subnames_fe)
                
                # Random effects design matrix
                X_list_re[[j]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
                if(ncol(X_list_re[[j]]) > 0) {
                    subnames_re <- paste0(par_name, ".", term_names[-(1:gam_setup$nsdf)])
                    names_re <- c(names_re, subnames_re)                    
                }
                
                # Smoothing matrix
                S_list[[j]] <- bdiag_check(gam_setup$S)
                
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
                        X_list_re = X_list_re, S_list = S_list,
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
        
        ###################
        ## Model fitting ##
        ###################
        #' @description TMB setup
        #'  
        #' @details This creates an attribute \code{tmb_obj}, which can be used to 
        #' evaluate the negative log-likelihood function.
        #' 
        #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
        #' @param map List passed to MakeADFun to fix parameters. (See TMB documentation.)
        setup = function(silent = TRUE, map = NULL) {
            # Number of time steps
            n <- nrow(self$data())
            
            # Create model matrices
            X_fe <- self$mats()$X_fe
            X_re <- self$mats()$X_re
            S <- self$mats()$S
            ncol_fe <- self$terms()$ncol_fe
            ncol_re <- self$terms()$ncol_re
            
            # Format initial parameters for TMB
            # (First fixed effects, then random effects)
            tmb_par <- list(coeff_fe = self$coeff_fe(),
                            log_lambda = 0,
                            log_decay = log(self$rho()),
                            coeff_re = 0)
            
            # Setup random effects
            random <- NULL
            if(is.null(S)) {
                # If there are no random effects, 
                # coeff_re and log_lambda are not estimated
                map <- c(map, list(coeff_re = factor(NA),
                                   log_lambda = factor(NA)))
                S <- as_sparse(matrix(0, 1, 1))
                ncol_re <- 0
                X_re <- as_sparse(rep(0, nrow(X_fe)))
            } else {
                # If there are random effects, 
                # set initial values for coeff_re and log_lambda
                random <- c(random, "coeff_re")
                tmb_par$coeff_re <- self$coeff_re()
                tmb_par$log_lambda <- log(self$lambda())
            }
            
            # TMB data object
            tmb_dat <- list(type = self$type(),
                            ID = self$data()$ID,
                            times = self$data()$time,
                            obs = as.matrix(self$obs()),
                            X_fe = as_sparse(X_fe),
                            X_re = as_sparse(X_re),
                            S = as_sparse(S),
                            ncol_re = ncol_re,
                            include_penalty = 1)
            
            # Model-specific data objects
            if(self$type() == "BM_t") {
                # Pass degrees of freedom for BM_t model
                tmb_dat$other_data <- self$other_data()$df
            } else if(self$type() == "BM_SSM" | self$type() == "OU_SSM") {
                # Number of dimensions
                n_dim <- ncol(self$obs())
                # Define initial state and covariance for Kalman filter
                # First index for each ID
                i0 <- c(1, which(self$data()$ID[-n] != self$data()$ID[-1]) + 1)
                # Initial state = first observation
                a0 <- as.matrix(self$obs()[i0,])
                tmb_dat$a0 <- a0
                # Initial state covariance
                if(is.null(self$other_data()$P0)) {
                    # Default if P0 not provided by user
                    tmb_dat$P0 <- diag(rep(10, n_dim))                    
                } else {
                    tmb_dat$P0 <- self$other_data()$P0
                }
                
                # Initialise model-specific parameter (measurement error SD)
                tmb_par <- c(log_sigma_obs =  0, tmb_par)
                
                # Check whether observation error is provided by user
                if(!is.null(self$other_data()$H)) {
                    tmb_dat$H_array <- self$other_data()$H
                    map <- c(map, list(log_sigma_obs = factor(NA)))
                } else {
                    tmb_dat$H_array <- array(0)
                }
            } else if(self$type() == "CTCRW") {
                # Number of dimensions
                n_dim <- ncol(self$obs())
                # Define initial state and covariance for Kalman filter
                # First index for each ID
                i0 <- c(1, which(self$data()$ID[-n] != self$data()$ID[-1]) + 1)
                # Initial state = (x1, 0, y1, 0, ...)
                a0 <- matrix(0, length(i0), 2*n_dim)
                for(i in 1:n_dim) {
                    a0[, 2*(i-1)+1] <- self$obs()[i0, i]
                }
                tmb_dat$a0 <- a0
                # Initial state covariance
                if(is.null(self$other_data()$P0)) {
                    # Default if P0 not provided by user
                    tmb_dat$P0 <- diag(rep(c(1, 10), n_dim))                    
                } else {
                    tmb_dat$P0 <- self$other_data()$P0
                }
                
                # Initialise model-specific parameter (measurement error SD)
                tmb_par <- c(log_sigma_obs =  0, tmb_par)
                
                # Check whether observation error is provided by user
                if(!is.null(self$other_data()$H)) {
                    tmb_dat$H_array <- self$other_data()$H
                    map <- c(map, list(log_sigma_obs = factor(NA)))
                } else {
                    tmb_dat$H_array <- array(0)
                }
            } else if(self$type() == "ESEAL_SSM") {
                # Define initial state and covariance for Kalman filter
                # Initial state = initial lipid mass
                tmb_dat$a0 <- cbind(1, rle(self$data()$dep_fat)$values)
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
            } else {
                # Unused for BM, OU, CIR...
                tmb_dat$other_data <- 0
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
            
            # Decaying response model
            if(self$type() %in% c("BM", "BM_t", "OU", "CIR")) {
                if(!is.null(self$other_data()$t_decay)) {
                    if(any(self$other_data()$col_decay > length(self$terms()$names_re_all))) {
                        stop(paste0("'col_decay' should be between 1 and ", 
                                    length(self$terms()$names_re_all)))
                    }
                    tmb_dat$t_decay <- self$other_data()$t_decay
                    tmb_dat$col_decay <- self$other_data()$col_decay
                    tmb_dat$ind_decay <- self$other_data()$ind_decay
                } else  {
                    tmb_dat$t_decay <- 0
                    tmb_dat$col_decay <- 0
                    tmb_dat$ind_decay <- 0
                    map <- c(map, list(log_decay = factor(NA)))
                }
            } else {
                # Remove log_decay if this is a model for which it isn't implemented
                tmb_par$log_decay <- NULL
            }
            
            # Create TMB object
            tmb_obj <- MakeADFun(data = tmb_dat, parameters = tmb_par, 
                                 DLL = "smoothSDE", silent = silent,
                                 map = map, random = random)
            
            # Negative log-likelihood function
            private$tmb_obj_ <- tmb_obj
            
            # Joint negative log-likelihood function excluding penalty
            # (used for conditional AIC)
            tmb_dat$include_penalty <- 0
            tmb_obj_joint <- MakeADFun(data = tmb_dat, parameters = tmb_par, 
                                       DLL = "smoothSDE", 
                                       map = map, silent = silent)
            private$tmb_obj_joint_ <- tmb_obj_joint
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
        #' @param map List passed to MakeADFun to fix parameters. See TMB documentation.
        fit = function(silent = TRUE, map = NULL) {
            # Print model formulation
            self$message()
            
            # Setup if necessary
            if(is.null(private$tmb_obj_)) {
                self$setup(silent = silent, map = map)
            }
            
            sys_time <- system.time({
                # Fit model
                private$out_ <- optim(par = private$tmb_obj_$par,
                                      fn = private$tmb_obj_$fn, 
                                      gr = private$tmb_obj_$gr)
                # private$out_ <- do.call(optim, private$tmb_obj_)
            })
            private$out_$systime <- sys_time
            # Get estimates and precision matrix for all parameters
            private$tmb_rep_ <- sdreport(private$tmb_obj_, 
                                         getJointPrecision = TRUE, 
                                         skip.delta.method = FALSE)
            
            # Save parameter estimates
            par_list <- as.list(private$tmb_rep_, "Estimate")
            self$update_coeff_fe(par_list$coeff_fe)
            if(length(self$terms()$ncol_re) > 0) {
                # Only save coeff_re and lambda is there are random effects
                self$update_coeff_re(par_list$coeff_re)
                self$update_lambda(exp(par_list$log_lambda))
            }
            
            if(!is.null(self$other_data()$t_decay)) {
                # Update decay rate parameter if decay model
                log_decay <- par_list$log_decay
                self$update_rho(exp(log_decay))
            }
        },
        
        ####################
        ## Get parameters ##
        ####################
        #' @description Get linear predictor for SDE parameters
        #' 
        #' @param new_data Optional data set of covariates. If \code{new_data},
        #' \code{X_fe} and \code{X_re} are not provided, then the observed 
        #' covariates are used.
        #' @param t Time points for which the parameters should be returned.
        #' If "all", returns parameters for all time steps (default).
        #' @param X_fe Optional design matrix for fixed effects, as returned
        #' by \code{make_mat}. If \code{new_data}, \code{X_fe} and \code{X_re} 
        #' are not provided, then the observed covariates are used.
        #' @param X_re Optional design matrix for random effects, as returned
        #' by \code{make_mat}. If \code{new_data}, \code{X_fe} and \code{X_re} 
        #' are not provided, then the observed covariates are used.
        #' @param coeff_fe Optional vector of fixed effect parameters
        #' @param coeff_re Optional vector of random effect parameters
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' @return Matrix of linear predictor 
        #' (X_fe %*% coeff_fe + X_re %*% coeff_re) 
        #' with one row for each time step and one column for each SDE parameter
        linear_predictor = function(new_data = NULL, t = "all",
                                    X_fe = NULL, X_re = NULL,
                                    coeff_fe = NULL, coeff_re = NULL,
                                    term = NULL) {
            # Get design matrices (X_fe/X_re) if not provided
            if(is.null(X_fe) | is.null(X_re)) {
                mats <- self$make_mat(new_data = new_data)
                if(is.null(X_fe)) {
                    X_fe <- mats$X_fe
                }
                if(is.null(X_re)) {
                    X_re <- mats$X_re
                }
            }
            
            # Use coeff_fe/coeff_re from model if not provided
            if(is.null(coeff_fe)) {
                coeff_fe <- self$coeff_fe()
            }
            if(is.null(coeff_re)) {
                coeff_re <- self$coeff_re()
            }
            
            # Only keep non-zero coefficients for relevant term
            if(!is.null(term)) {
                coeff_fe_term <- rep(0, length(coeff_fe))
                coeff_re_term <- rep(0, length(coeff_re))
                term_ind <- term_indices(names_fe = self$terms()$names_fe, 
                                         names_re = self$terms()$names_re_all, 
                                         term = term)
                coeff_fe_term[term_ind$fe] <- coeff_fe[term_ind$fe]
                coeff_re_term[term_ind$re] <- coeff_re[term_ind$re]
                coeff_fe <- coeff_fe_term
                coeff_re <- coeff_re_term
            }
            
            # Get linear predictor and format into matrix
            lp <- X_fe %*% coeff_fe + X_re %*% coeff_re
            lp_mat <- matrix(lp, ncol = length(self$formulas()))
            colnames(lp_mat) <- names(self$formulas())
            
            # Keep rows of lp_mat given in 't'
            if(length(t) == 1) { if(t == "all") {
                t <- 1:nrow(lp_mat)
            }}
            if(any(t < 1 | t > nrow(lp_mat))) {
                stop("Elements of 't' should be between 1 and", nrow(lp_mat))
            }
            lp_mat <- lp_mat[t,, drop = FALSE]
            
            return(lp_mat)
        },
        
        #' @description Get SDE parameters
        #' 
        #' @param t Time points for which the parameters should be returned.
        #' If "all", returns parameters for all time steps. Default: 1.
        #' @param new_data Optional data set of covariates. If \code{new_data},
        #' \code{X_fe} and \code{X_re} are not provided, then the observed 
        #' covariates are used.
        #' @param X_fe Optional design matrix for fixed effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param X_re Optional design matrix for random effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param coeff_fe Optional vector of fixed effect parameters
        #' @param coeff_re Optional vector of random effect parameters
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' @return Matrix with one row for each time point in t, and one
        #' column for each SDE parameter
        par = function(t = NULL, new_data = NULL, 
                       X_fe = NULL, X_re = NULL,
                       coeff_fe = NULL, coeff_re = NULL, 
                       resp = TRUE, term = NULL) {
            # Default t = 1, unless new data are provided (then t = all)
            if(is.null(t)) {
                if(!is.null(new_data) | !is.null(X_fe) | !is.null(X_re)) {
                    t <- "all"
                } else {
                    t <- 1   
                }
            }
            
            # Get linear predictor
            lp_mat <- self$linear_predictor(new_data = new_data, t = t, 
                                            X_fe = X_fe, X_re = X_re,
                                            coeff_fe = coeff_fe, 
                                            coeff_re = coeff_re, 
                                            term = term)
            
            # Apply inverse link to get parameters on response scale
            if(resp) {
                par_mat <- matrix(sapply(1:ncol(lp_mat), function(i) {
                    self$invlink()[[i]](lp_mat[,i])   
                }), ncol = ncol(lp_mat))
                colnames(par_mat) <- names(self$invlink())
            } else {
                par_mat <- lp_mat
            }
            return(par_mat)
        },
        
        ################################
        ## Uncertainty quantification ##
        ################################
        #' @description Posterior draws (coefficients)
        #' 
        #' @param n_post Number of posterior draws
        #' 
        #' @return Matrix with one column for each coefficient and one row for each
        #' posterior draw
        post_coeff = function(n_post) {
            # Number of SDE parameters
            n_par <- length(self$formulas())
            # TMB report
            rep <- self$tmb_rep()
            
            # Joint covariance matrix
            if(!is.null(rep$jointPrecision)) {
                jointCov <- prec_to_cov(rep$jointPrecision)
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
            post_coeff <- matrix(rmvn(n = n_post, mu = par_all, V = jointCov), 
                                 ncol = ncol(jointCov))
            
            # Split matrix into list of matrices
            names_coeff <- colnames(jointCov)
            post_list <- lapply(unique(names_coeff), function(name) 
                post_coeff[, which(names_coeff == name), drop = FALSE])
            names(post_list) <- unique(names_coeff)
            
            # Add empty vector for coeff_re if no random effects
            if(!"coeff_re" %in% unique(names_coeff)) {
                post_list$coeff_re <- matrix(NA, nrow = n_post, ncol = 0)
            }
            
            # Deal with fixed SDE parameters
            ind_estpar <- which(!names(self$formulas()) %in% self$fixpar())
            # Indices of estimated coefficients in coeff_fe
            fe_cols <- rep(1:n_par, self$terms()$ncol_fe)
            ind_est_fe <- which(fe_cols %in% ind_estpar)
            
            # In post_fe, set columns for fixed parameters to fixed value,
            # and use posterior draws for non-fixed parameters
            post_fe <- matrix(rep(self$coeff_fe(), each = n_post),
                              nrow = n_post, ncol = sum(self$terms()$ncol_fe))
            post_fe[,ind_est_fe] <- post_list$coeff_fe
            post_list$coeff_fe <- post_fe
            
            # Set column names
            colnames(post_list$coeff_fe) <- self$terms()$names_fe
            colnames(post_list$coeff_re) <- self$terms()$names_re_all
            
            return(post_list)
        },
        
        #' @description Posterior draws of SDE parameters (for uncertainty 
        #' quantification)
        #' 
        #' @param X_fe Design matrix (fixed effects)
        #' @param X_re Design matrix (random effects)
        #' @param n_post Number of posterior draws (default: 100)
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' @return Array with one row for each time step, one column for
        #' each SDE parameter, and one layer for each posterior draw
        post_par = function(X_fe, X_re, n_post = 100, resp = TRUE, term = NULL) {
            # Number of SDE parameters
            n_par <- length(self$formulas())
            # Number of time steps
            n <- nrow(X_fe)/n_par
            
            # Generate posterior draws of coeff_fe and coeff_re
            post_coeff <- self$post_coeff(n_post = n_post)
            
            # Get corresponding SDE parameters
            par_array <- array(sapply(1:n_post, function(i) {
                self$par(t = "all", 
                         X_fe = X_fe, X_re = X_re, 
                         coeff_fe = post_coeff$coeff_fe[i,], 
                         coeff_re = post_coeff$coeff_re[i,],
                         resp = resp, 
                         term = term)  
            }), dim = c(n, n_par, n_post))
            dimnames(par_array)[[2]] <- names(self$invlink())
            
            return(par_array)
        },
        
        #' @description Pointwise confidence intervals for SDE parameters
        #'
        #' @param t Time points for which the parameters should be returned.
        #' If "all", returns parameters for all time steps. Defaults to 1 if new
        #' data not provided, or "all" if new data provided.
        #' @param new_data Optional data frame containing covariate values 
        #' for which the CIs should be computed
        #' @param X_fe Optional design matrix for fixed effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param X_re Optional design matrix for random effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param level Confidence level (default: 0.95 for 95\% confidence 
        #' intervals)
        #' @param n_post Number of posterior samples from which the confidence
        #' intervals are calculated. Larger values will reduce approximation
        #' error, but increase computation time. Defaults to 1000.
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' This method generates pointwise confidence intervals 
        #' by simulation. That is, it generates \code{n_post} posterior samples 
        #' of the estimated parameters from a multivariate normal distribution,
        #' where the mean is the vector of estimates and the covariance matrix 
        #' is provided by TMB (using \code{\link{post_par}}). Then, the SDE 
        #' parameters are derived for each set of posterior parameter values, 
        #' and pointwise confidence intervals are obtained as quantiles of the 
        #' posterior simulated SDE parameters.
        #' 
        #' @return List with elements:
        #' \itemize{
        #'   \item{\code{low}}{Matrix of lower bounds of confidence intervals.}
        #'   \item{\code{upp}}{Matrix of upper bounds of confidence intervals.}
        #' }
        CI_pointwise = function(t = NULL, new_data = NULL, 
                                X_fe = NULL, X_re = NULL, 
                                level = 0.95, n_post = 1e3, 
                                resp = TRUE, term = NULL) {
            # Default t = 1, unless new data are provided (then t = all)
            if(is.null(t)) {
                if(!is.null(new_data) | !is.null(X_fe) | !is.null(X_re)) {
                    t <- "all"
                } else {
                    t <- 1   
                }
            }
            
            if(is.null(X_fe) | is.null(X_re)) {
                if(is.null(new_data)) {
                    new_data <- self$data()
                }
                if(t != "all") {
                    new_data <- new_data[t,]
                }
                mats <- self$make_mat(new_data = new_data)
                X_fe <- mats$X_fe
                X_re <- mats$X_re
            }
            
            # Posterior samples of SDE parameters
            post_par <- self$post_par(X_fe = X_fe, X_re = X_re, 
                                      n_post = n_post, resp = resp,
                                      term = term)
            
            # Get confidence intervals as quantiles of posterior samples
            alpha <- (1 - level)/2
            CI <- apply(post_par, c(1, 2), quantile, probs = c(alpha, 1 - alpha))
            
            # Array with one row per parameter, columns for lower/upper bounds,
            # and one layer for each time step
            CI <- aperm(CI, c(3, 1, 2))
            dimnames(CI)[[2]] <- c("low", "upp")
            
            return(CI)
        },
        
        #' @description Simultaneous confidence intervals for SDE parameters
        #' 
        #' @param t Time points for which the parameters should be returned.
        #' If "all", returns parameters for all time steps. Defaults to 1 if new
        #' data not provided, or "all" if new data provided.
        #' @param new_data Optional data frame containing covariate values 
        #' for which the CIs should be computed
        #' @param X_fe Optional design matrix for fixed effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param X_re Optional design matrix for random effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param level Confidence level (default: 0.95 for 95\% confidence 
        #' intervals)
        #' @param n_post Number of posterior samples from which the confidence
        #' intervals are calculated. Larger values will reduce approximation
        #' error, but increase computation time. Defaults to 1000.
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' This method closely follows the approach suggested by Gavin Simpson at
        #' fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/,
        #' itself based on Section 6.5 of Ruppert et al. (2003).
        #' 
        #' @return List with elements:
        #' \itemize{
        #'   \item{\code{low}}{Matrix of lower bounds of confidence intervals.}
        #'   \item{\code{upp}}{Matrix of upper bounds of confidence intervals.}
        #' }
        CI_simultaneous = function(t = NULL, new_data = NULL, 
                                   X_fe = NULL, X_re = NULL, 
                                   level = 0.95, n_post = 1000, 
                                   resp = TRUE, term = NULL) {
            # Default t = 1, unless new data are provided (then t = all)
            if(is.null(t)) {
                if(!is.null(new_data) | !is.null(X_fe) | !is.null(X_re)) {
                    t <- "all"
                } else {
                    t <- 1   
                }
            }
            
            if(is.null(X_fe) | is.null(X_re)) {
                if(is.null(new_data)) {
                    new_data <- self$data()
                }
                if(t != "all") {
                    new_data <- new_data[t,]
                }
                mats <- self$make_mat(new_data = new_data)
                X_fe <- mats$X_fe
                X_re <- mats$X_re
            }
            
            # Number of parameters of this SDE model
            n_par <- length(self$formulas())
            # Number of time steps
            n <- nrow(X_fe)/n_par
            
            # Get SE of parameters on linear predictor scale
            par_linpred <- self$par(t = "all", X_fe = X_fe, X_re = X_re, 
                                    resp = FALSE, term = term)
            CIpw_linpred <- self$CI_pointwise(X_fe = X_fe, X_re = X_re, 
                                              level = level, n_post = n_post, 
                                              resp = FALSE, term = term)
            se_linpred <- (par_linpred - t(CIpw_linpred[,"low",]))/qnorm((1 + level)/2)
            se_linpred_vec <- as.vector(se_linpred)
            
            # Posterior samples of (estimate - true par) ~ N(0, V)
            coeff_post <- self$post_coeff(n_post = n_post)
            diff_fe <- t(sweep(x = coeff_post$coeff_fe, MARGIN = 2, 
                               STATS = self$coeff_fe(), FUN = "-"))
            diff_re <- t(sweep(x = coeff_post$coeff_re, MARGIN = 2, 
                               STATS = self$coeff_re(), FUN = "-"))
            
            # sim_dev <- sapply(1:n_post, function(i) {
            #     lp <- self$linear_predictor(
            #         new_data = new_data, t = "all",
            #         coeff_fe = diff_fe[,i], coeff_re = diff_re[,i], 
            #         term = term)
            #     return(as.vector(lp))
            # })
            
            # The following is redundant with linear_predictor because, at the
            # moment, calling linear_predictor to each column of diff_fe/re is
            # too computational. In the future, linear_predictor should be able
            # to work with matrices provided for coeff_fe/re.
            
            # Only keep non-zero coefficients for relevant term
            if(!is.null(term)) {
                diff_fe_term <- matrix(0, nrow = nrow(diff_fe), ncol = ncol(diff_fe))
                diff_re_term <- matrix(0, nrow = nrow(diff_re), ncol = ncol(diff_re))
                term_ind <- term_indices(names_fe = self$terms()$names_fe, 
                                         names_re = self$terms()$names_re_all, 
                                         term = term)
                diff_fe_term[term_ind$fe,] <- diff_fe[term_ind$fe,]
                diff_re_term[term_ind$re,] <- diff_re[term_ind$re,]
                diff_fe <- diff_fe_term
                diff_re <- diff_re_term
            }
            
            # Deviations between fitted and true parameters
            sim_dev <- X_fe %*% diff_fe + X_re %*% diff_re
            
            # Absolute values of the standardized deviations from the true model
            abs_dev <- abs(sweep(sim_dev, 1, se_linpred_vec, FUN = "/"))
            abs_dev_array <- array(abs_dev, dim = c(n, n_par, n_post))
            
            # Replace 0/0 = NaN by zero (case where there is no deviation)
            abs_dev_array[which(is.nan(abs_dev_array))] <- 0
            
            # Take maximum
            max_abs_dev <- apply(abs_dev_array, c(2, 3), max)
            
            # Critical value (na.rm for case where SE = 0 and abs_dev = NaN)
            crit <- apply(max_abs_dev, 1, quantile, prob = level, na.rm = TRUE)
            crit[which(is.na(crit))] <- 0
            
            # Lower and upper bounds of simultaneous CIs
            dimnames <- list(NULL, c("low", "upp"), names(self$formulas()))
            CIs <- array(sapply(1:n_par, function(i) {
                par <- names(self$formulas())[i]
                invlink <- ifelse(resp, yes = self$invlink()[[par]], no = identity)
                low <- invlink(par_linpred[,par] - (crit[i] * se_linpred[,par]))
                upp <- invlink(par_linpred[,par] + (crit[i] * se_linpred[,par]))
                return(matrix(c(low, upp), ncol = 2))
            }), dim = c(n, 2, n_par),  dimnames = dimnames)
            
            CIs <- aperm(CIs, c(3, 2, 1))
            return(CIs)
        },
        
        ####################
        ## Model checking ##
        ####################
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
            par <- self$par(t = "all")
            
            # Response variable(s)
            Z <- as.matrix(data[, self$response(), drop = FALSE])
            
            # Compute mean and sd of normal transition density
            if(self$type() == "BM") {
                mean <- Z[-end_ind,] + par[-end_ind, "mu"] * dtimes
                sd <- par[-end_ind, "sigma"] * sqrt(dtimes)
            } else if (self$type() == "BM_t") {
                df <- self$other_data()$df # degrees of freedom
                mean <- Z[-end_ind,] + par[-end_ind, "mu"] * dtimes
                sd <- par[-end_ind, "sigma"] * sqrt(dtimes)
                sd <- sd / sqrt(df / (df-2)) # this is actually the scale parameter, not the SD any more
            } else if(self$type() == "OU") {
                # Need to account for case where mu has several dimensions
                mu <- par[-end_ind, which(substr(colnames(par), 1, 2) == "mu")]
                tau <- par[-end_ind, "tau"]
                kappa <- par[-end_ind, "kappa"]
                mean <-  mu + exp(- dtimes/tau) * (Z[-end_ind,] - mu)
                sd <- sqrt(kappa * (1 - exp(-2 * dtimes / tau)))
            } else {
                stop(paste("Residuals not implemented for model", self$type()))
            }
            
            # Residuals ~ N(0, 1) under assumptions of model and discretization
            res <- matrix(NA, nrow = n, ncol = ncol(Z))
            res[-end_ind,] <- (Z[-start_ind,] - mean) / sd
            return(res)
        },
        
        #' @description Posterior predictive checks
        #' 
        #' @param check_fn Goodness-of-fit function which accepts "data" as input
        #' and returns a statistic or a vector of statistics, to be compared
        #' between observed data and simulations. 
        #' @param n_sims Number of simulations to perform 
        #' @param silent Logical. If FALSE, simulation progress is shown. 
        #' (Default: FALSE)
        #' 
        #' @details This applies \code{check_fn} to the observed data (returned by 
        #' \code{data()} method) to obtain observed statistics. It then repeatedly
        #' simulates a realisation from the fitted SDE (based on observed covariates),
        #' and applies \code{check_fn} to the simulated data. The simulations use
        #' the \code{posterior = TRUE} option from the code{simulate} method, i.e., 
        #' parameters of model used for simulation are generated from posterior
        #' distribution.
        #'
        #' @return List with elements:
        #' \itemize{
        #'   \item{obs_stat}{Vector of values of goodness-of-fit statistics for the
        #'   observed data}
        #'   \item{stats}{Matrix of values of goodness-of-fit statistics for the
        #'   simulated data sets (one row for each statistic, and one column for each
        #'   simulation)}
        #'   \item{plot}{ggplot object showing, for each statistic returned by 
        #'   check_fn, a histogram of simulated values and a vertical line for the
        #'   observed value. An observed value in the tails of the simulated
        #'   distribution suggests lack of fit.}
        #' }
        check_post = function(check_fn, n_sims = 100, silent = FALSE) {
            # Evaluate statistics for observed data
            obs_stat <- check_fn(self$data())
            
            # Simulate from model and evaluate statistics for simulated data
            stats <- matrix(0, nrow = length(obs_stat), ncol = n_sims)
            for (sim in 1:n_sims) {
                if (!silent) {
                    cat(paste0("Simulation ", sim, "/", n_sims, "\r"))
                }
                
                # Simulate new data
                new_data <- self$simulate(data = self$data(), posterior = TRUE) 
                # Compute statistics
                stats[,sim] <- check_fn(new_data)
            }
            
            # Get names of statistics
            if(is.null(names(obs_stat))) {
                stat_names <- paste("statistic", 1:length(obs_stat))
                names(obs_stat) <- stat_names
            } else {
                stat_names <- names(obs_stat)
            }
            rownames(stats) <- stat_names
            
            # Data frames for ggplot (observed and simulated values)
            df_obs <- as.data.frame.table(as.matrix(obs_stat))
            colnames(df_obs) <- c("stat", "sim", "val")
            df <- as.data.frame.table(stats)
            colnames(df) <- c("stat", "sim", "val")
            
            # Create plot
            p <- ggplot(df, aes(val)) + 
                geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                               bg = "lightgrey", na.rm = TRUE) +
                facet_wrap("stat", scales = "free") + 
                geom_vline(aes(xintercept = val), data = df_obs) +
                xlab("statistic") + ggtitle("Vertical line is observed value") +
                theme_light() +
                theme(strip.background = element_blank(),
                      strip.text = element_text(colour = "black"))
            if (!silent) {
                plot(p)
            }
            
            return(list(obs_stat = obs_stat, stats = stats, plot = p))
        }, 
        
        #' Conditional Akaike Information Criterion
        #' 
        #' The conditional AIC is for example defined by 
        #' Wood (2017), as AIC = - 2L + 2k where L is the
        #' maximum joint log-likelihood (of fixed and random
        #' effects), and k is the number of effective degrees
        #' of freedom of the model (accounting for flexibility
        #' in non-parametric terms implied by smoothing)
        #' 
        #' @return Conditional AIC
        AIC_conditional = function() {
            # Effective degrees of freedom
            edf <- self$edf_conditional()
            
            # Joint likelihood
            par_all <- c(self$tmb_rep()$par.fixed, self$tmb_rep()$par.random)
            llk <- - private$tmb_obj_joint_$fn(par_all)
            
            aic <- - 2 * llk + 2 * edf
            return(aic)
        },
        
        #' Marginal Akaike Information Criterion
        #' 
        #' The marginal AIC is for example defined by 
        #' Wood (2017), as AIC = - 2L + 2k where L is the
        #' maximum marginal log-likelihood (of fixed 
        #' effects), and k is the number of degrees
        #' of freedom of the fixed effect component of
        #' the model
        #' 
        #' @return Marginal AIC
        AIC_marginal = function() {
            # Fixed effect DF
            edf <- length(self$out()$par) - length(self$lambda())
            
            # Marginal likelihood
            llk <- - self$out()$value
            
            aic <- - 2 * llk + 2 * edf
            return(aic)
        },
        
        #' @description Effective degrees of freedom
        #'
        #' @return Number of effective degrees of freedom
        #' (accounting for flexibility in non-parametric 
        #' terms implied by smoothing)
        edf_conditional = function() {
            # Degrees of freedom for fixed effects
            edf <- length(self$out()$par) - length(self$lambda())
            
            if(!is.null(self$tmb_rep()$jointPrecision)) {
                # Get Hessian of joint neg log-likelihood
                par_all <- c(self$tmb_rep()$par.fixed, self$tmb_rep()$par.random)
                H <- self$tmb_obj_joint()$he(par_all)
                
                # Joint covariance matrix
                Q <- self$tmb_rep()$jointPrecision
                V <- prec_to_cov(Q)
                
                # Subset to random effect components
                ind_re <- which(colnames(Q) == "coeff_re")
                V_re <- V[ind_re, ind_re]
                H_re <- H[ind_re, ind_re]
                
                # Random effect EDF
                edf <- edf + sum(diag(H_re %*% V_re))
            }
            
            return(edf)
        },
        
        ################
        ## Simulation ##
        ################
        #' @description Simulate from SDE model
        #' 
        #' @param data Data frame for input data. Should have at least one column 'time' for
        #' times of observations, and columns for covariates if necessary.
        #' @param z0 Optional value for first observation of simulated time series.
        #' Default: 0.
        #' @param posterior Logical. If TRUE, the parameters used for the simulation
        #' are drawn from their posterior distribution using \code{SDE$post_coeff}, 
        #' therefore accounting for uncertainty.
        #' 
        #' @return Input data frame with extra column for simulated time series
        simulate = function(data, z0 = 0, posterior = FALSE) {
            # Check that data includes times of observations
            if(is.null(data$time)) {
                stop("'data' should have a column named 'time'")
            }
            if(is.null(data$ID)) {
                data$ID <- factor(1)
            }
            
            # Create SDE parameters
            if(posterior) {
                # Use posterior draw for coeff_fe/re
                coeff <- self$post_coeff(n_post = 1)
                par <- self$par(new_data = data, 
                                coeff_fe = coeff$coeff_fe[1,],
                                coeff_re = coeff$coeff_re[1,])
            } else {
                # Use point estimates of coeff_fe/re
                par <- self$par(new_data = data)
            }
            
            # Loop over dimensions
            n_dim <- length(self$response())
            if(length(z0) < n_dim) {
                z0 <- rep(z0, n_dim)
            }
            for(d in 1:n_dim) {
                # Initialize vector of simulated observations
                obs <- rep(NA, nrow(data))
                
                # Loop over IDs
                for(id in seq_along(unique(data$ID))) {
                    # Get relevant rows of data
                    ind <- which(data$ID == unique(data$ID)[id])
                    dtimes <- diff(data$time[ind])
                    sub_n <- length(ind)
                    sub_obs <- rep(z0[d], sub_n)
                    sub_par <- par[ind,]
                    
                    if(self$type() == "BM") {
                        # If BM, generate all increments directly
                        mean <- sub_par[-sub_n, d] * dtimes
                        sd <- sub_par[-sub_n, n_dim + 1] * sqrt(dtimes)
                        sub_obs <- cumsum(c(z0[d], rnorm(sub_n - 1, mean = mean, sd = sd)))
                    } else if(self$type() == "OU") {
                        # If OU, loop over observation times
                        for(i in 2:sub_n) {
                            # Generate observation from OU transition density
                            mean <- exp(- dtimes[i-1] / sub_par[i-1, n_dim + 1]) * sub_obs[i-1] +
                                (1 - exp(-dtimes[i-1] / sub_par[i-1, n_dim + 1])) * sub_par[i-1, d]
                            sd <- sqrt(sub_par[i-1,  n_dim + 2] * 
                                           (1 - exp(-2 * dtimes[i-1] / sub_par[i-1, n_dim + 1])))
                            sub_obs[i] <- rnorm(n = 1, mean = mean, sd = sd)
                        }
                    } else if(self$type() == "CTCRW") {
                        # Data has one column for velocity and one for position
                        sub_dat <- matrix(rep(c(0, z0[d]), each = sub_n), nrow = sub_n,
                                          dimnames = list(NULL, c("v", "z")))
                        
                        # Unpack parameters
                        mu <- sub_par[, d]
                        tau <- sub_par[, n_dim + 1]
                        nu <- sub_par[, n_dim + 2]
                        beta <- 1/tau
                        sigma <- 2 * nu / sqrt(tau * pi)
                        
                        # Loop over time steps
                        mean <- rep(NA, 2)
                        for(i in 2:sub_n) {
                            # Mean of next state vector (V, Z)
                            p <- exp(-beta[i-1] * dtimes[i-1])
                            mean[1] <- p * sub_dat[i-1, "v"] + (1 - p) * mu[i-1]
                            mean[2] <- sub_dat[i-1, "z"] + mu[i-1] * dtimes[i-1] +
                                (sub_dat[i-1, "v"] - mu[i-1]) / beta[i-1] * (1 - p)
                            
                            # Covariance of next state vector
                            V <- CTCRW_cov(beta = beta[i-1], sigma = sigma[i-1], 
                                           dt = dtimes[i-1])
                            
                            sub_dat[i,] <- rmvn(1, mu = mean, V = V)
                        }
                        
                        # Only return location (and not velocity)
                        sub_obs <- sub_dat[,"z"]
                    } else if(self$type() == "CIR") {
                        # Unpack parameters
                        mu <- sub_par[, d]
                        beta <- sub_par[, n_dim + 1]
                        sigma <- sub_par[, n_dim + 2]
                        
                        # Loop over time steps
                        sub_obs <- rep(z0[d], sub_n)
                        for(i in 2:n) {
                            c <- 2 * beta[i-1] / 
                                ((1 - exp(-beta[i-1] * dtimes[i-1])) * sigma[i-1]^2)
                            df <- 4 * beta[i-1] * mu[i-1] / sigma[i-1]^2
                            ncp <- 2 * c * sub_obs[i-1] * exp(- beta * dtimes[i-1])
                            Y <- rchisq(n = 1, df = df, ncp = ncp)
                            sub_obs[i] <- Y/(2 * c)
                        }
                    } else {
                        stop(paste("Simulation not implemented yet for", self$type(), "model."))
                    }
                    
                    # Update observation vector
                    obs[ind] <- sub_obs
                }
                
                # Add simulated variable to data frame
                data[[self$response()[d]]] <- obs
            }
            
            return(data)
        },
        
        ######################
        ## Plotting methods ##
        ######################
        #' @description Plot observation parameters
        #' 
        #' @param var Name of covariate as a function of which the parameters
        #' should be plotted
        #' @param par_names Optional vector for the names of SDE parameters to plot.
        #' If not specified, all parameters are plotted.
        #' @param covs Optional data frame with a single row and one column
        #' for each covariate, giving the values that should be used. If this is
        #' not specified, the mean value is used for numeric variables, and the
        #' first level for factor variables.
        #' @param n_post Number of posterior draws to plot. Default: 100.
        #' @param show_CI Should confidence bands be plotted rather than posterior
        #' draws? Can takes values 'none' (default; no confidence bands), 
        #' 'pointwise' (show pointwise confidence bands obtained with
        #' \code{\link{CI_pointwise}}), or 'simultaneous' (show simultaneous 
        #' confidence bands obtained with \code{\link{CI_simultaneous}})
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param term Name of model term as character string, e.g., "time", 
        #' or "s(time)". Use \code{$coeff_fe()} and \code{$coeff_re()} methods
        #' to find names of model terms. This uses fairly naive substring 
        #' matching, and may not work if one covariate's name is a 
        #' substring of another one.
        #' 
        #' @return A ggplot object
        plot_par = function(var, par_names = NULL, covs = NULL, n_post = 100,
                            show_CI = "none", resp = TRUE, term = NULL) {
            if(missing(var)) {
                var_names <- unique(rapply(self$formulas(), all.vars))
                error_message <- paste0("'var' should be one of: '", 
                                        paste0(var_names, collapse = "', '"), "'")
                stop(error_message)
            }
            
            # Create design matrices
            mats <- self$make_mat_grid(var = var, covs = covs)
            par <- self$par(t = "all", X_fe = mats$X_fe, X_re = mats$X_re, 
                            resp = resp, term = term)
            
            # Data frame for posterior draws
            post_df <- NULL
            if(n_post > 0 & show_CI == "none") {
                post <- self$post_par(X_fe = mats$X_fe, X_re = mats$X_re, 
                                      n_post = n_post, resp = resp, term = term)
                post_df <- as.data.frame.table(post)
                colnames(post_df) <- c("var", "par", "stratum", "val")
                post_df$mle = "no"                   
            } else if (show_CI != "none") {
                CI_fn <- ifelse(show_CI == "pointwise", 
                                yes = self$CI_pointwise, 
                                no = self$CI_simultaneous)
                CI <- CI_fn(X_fe = mats$X_fe, X_re = mats$X_re,
                            n_post = n_post, level = 0.95,
                            resp = resp, term = term)
                CI_df <- data.frame(var = mats$new_data[,var],
                                    par = rep(names(self$formulas()), each = dim(CI)[3]),
                                    stratum = "CI",
                                    low = as.vector(t(CI[,"low",])),
                                    upp = as.vector(t(CI[,"upp",])))
            }
            
            # Data frame for MLE
            mle_df <- as.data.frame.table(par)
            colnames(mle_df) <- c("var", "par", "val")
            mle_df$stratum <- "mle"
            mle_df$mle <- "yes"
            
            # Full data frame
            df <- rbind(mle_df, post_df)
            df$var <- mats$new_data[, var]
            
            # If 'par_names' is specified, subset to relevant parameters
            if(!is.null(par_names)) {
                if(!all(par_names %in% unique(df$par))) {
                    error_message <- paste0("Check that elements of 'par_names' are in: '", 
                                            paste0(unique(df$par), collapse = "', '"), "'")
                    stop(error_message)
                }
                df <- subset(df, par %in% par_names)
                if(show_CI != "none") {
                    CI_df <- subset(CI_df, par %in% par_names)
                }
            }
            
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
            p0 <- ggplot(df, aes(x = var, group = stratum)) + 
                scale_colour_manual(values = pal, guide = "none") +
                facet_wrap(c("par"), scales = "free_y",
                           strip.position = "left",
                           labeller = label_bquote(.(as.character(par)))) +
                theme_light() + xlab(var) + ylab(NULL) + ggtitle(plot_txt) +
                theme(strip.background = element_blank(),
                      strip.placement = "outside", 
                      strip.text = element_text(colour = "black"))
            
            if(show_CI != "none") {
                p0 <- p0 + 
                    geom_ribbon(aes(ymin = low, ymax = upp), data = CI_df, 
                                fill = rgb(0.2, 0.5, 0.8, 0.3))
            }
            
            if(is.factor(df$var)) {
                p <- p0 +
                    geom_point(aes(y = val, size = mle, col = mle)) +
                    scale_size_manual(values = c(0.1, 1), guide = "none") +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
            } else {
                p <- p0 +
                    geom_line(aes(y = val, col = mle))
            }
            
            return(p)
        },
        
        ###################
        ## Miscellaneous ##
        ###################
        #' @description Indices of fixed coefficients in coeff_fe
        ind_fixcoeff = function() {
            # Number of SDE parameters
            n_par <- length(self$formulas())
            
            # Number of columns of X_fe for each SDE parameter
            ncol_fe <- self$terms()$ncol_fe
            
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
        
        #' @description Print equation for this model
        eqn = function() {
            switch (self$type(),
                    "BM" = "    dZ(t) = mu dt + sigma dW(t)",
                    "BM_SSM" = paste0("    dY(t) = mu dt + sigma dW(t)\n",
                                      "    Z(i) ~ N(Y(i), sigma_obs^2)"),
                    "BM_t" = "    Brownian motion with t-distributed noise",
                    "OU" = paste0("    dZ(t) = beta (mu - Z(t)) dt + sigma dW(t)\n",
                                  "Parameterised in terms of:\n",
                                  "* tau = 1/beta\n",
                                  "* kappa = sigma^2/(2*beta)"),
                    "OU_SSM" = paste0("    dZ(t) = beta (mu - Z(t)) dt + sigma dW(t)\n",
                                      "    Z(i) ~ N(Y(i), sigma_obs^2)\n",
                                      "Parameterised in terms of:\n",
                                      "* tau = 1/beta\n",
                                      "* kappa = sigma^2/(2*beta)"),
                    "CIR" = "    dZ(t) = beta (mu - Z(t)) dt + sigma sqrt(Z(t)) dW(t)",
                    "CTCRW" = paste0("    dV(t) = beta (mu - V(t)) dt + sigma dW(t)\n", 
                                     "    dZ(t) = V(t) dt\n",
                                     "Parameterised in terms of:\n",
                                     "* tau = 1/beta\n",
                                     "* nu = sqrt(pi/beta)*sigma/2"),
                    "ESEAL_SSM" = paste0("    dL(t) = mu dt + sigma dW(t)\n", 
                                         "    Z(i) ~ N(a1 + a2 L(i)/R(i), tau^2/h(i))"))
            
        },
        
        #' @description Print SDE and parameter formulas
        message = function() {
            message("#######################")
            message("### smoothSDE model ###")
            message("#######################")
            
            # Print SDE
            eqn <- self$eqn()
            message("> SDE for ", self$type(), " model:")
            message(eqn, "\n")
            
            # Print parameter formulas
            message("> Formulas for model parameters:")
            f <- self$formulas()
            for(i in 1:length(f)) {
                if(names(f)[i] %in% self$fixpar()) {
                    this_form <- "fixed"
                } else {
                    this_form <- as.character(f[[i]])[2]
                    this_form <- gsub("\\+", "+\n\t", this_form)
                }
                message("* ", names(f)[i], " ~ ", this_form)
            }
            message()
        },
        
        #' @description Print parameter values for t = 1
        print_par = function() {
            if(is.null(private$out_)) {
                message("> Initial SDE parameters (t = 1):")
            } else {
                message("> Estimated SDE parameters (t = 1):")
                CI <- round(self$CI_pointwise(t = 1), 3)
            }
            par <- self$par(t = 1)
            
            for(i in 1:length(par)) {
                msg <- paste0("* ", colnames(par)[i], " = ",  round(par[i], 3))
                if(!is.null(private$out_)) {
                    msg <- paste0(msg, "\t (", CI[i, 1,], ", ", CI[i, 2,], ")")
                }
                message(msg)
            }
        },
        
        #' @description Print SDE object
        print = function() {
            self$message()
            self$print_par()
        },
        
        #' @description Print stationary distribution
        stationary = function() {
            if(is.null(private$out_)) {
                msg <- "Based on initial SDE parameters (t = 1), "
            } else {
                msg <- "Based on estimated SDE parameters (t = 1), "
                CI <- round(self$CI_pointwise(t = 1), 3)
            }
            par <- round(self$par(t = 1), 3)
            msg <- paste0(msg, "the stationary distribution of this ", 
                          self$type(), " process is ")
            
            if(self$type() == "OU" | self$type() == "OU_SSM") {
                # OU: Z_t ~ N(mu, kappa)
                msg <- paste0(msg, "normal with parameters:\n",
                              "\t* mean = ", par[,"mu"], " \t(", 
                              CI["mu", 1,], ", ", CI["mu", 2,], ")\n",
                              "\t* variance = ", par[,"kappa"], " \t(", 
                              CI["kappa", 1,], ", ", CI["kappa", 2,], ")")
            } else if(self$type() == "CIR") {
                # Shape and rate of stationary distribution
                mean <- par[,"mu"]
                var <- par[,"mu"] * par[,"sigma"]^2 / (2 * par[,"beta"])
                
                # Generate posterior draws from SDE parameters to get CIs
                mats <- self$make_mat(new_data = self$data()[1,])
                post <- self$post_par(X_fe = mats$X_fe, X_re = mats$X_re, n_post = 1e3)
                post_mean <- post[,"mu",]
                post_var <- post[,"mu",] * post[,"sigma",]^2 / (2 * post[,"beta",])
                CI_mean <- round(quantile(post_mean, c(0.025, 0.975)), 3)
                CI_var <- round(quantile(post_var, c(0.025, 0.975)), 3)
                
                msg <- paste0(msg, "gamma with parameters:\n",
                              "\t* mean = ", round(mean, 3), " \t(",
                              CI_mean[1], ", ", CI_mean[2], ")\n",
                              "\t* variance = ", round(var, 3), " \t(",
                              CI_var[1], ", ", CI_var[2], ")")
            }
            
            msg <- paste0(msg, "\n(Note: this is *not* the stationary distribution ",
                          "if the parameters are time-varying)")
            message(msg)
        }
    ),
    
    private = list(
        formulas_ = NULL,
        data_ = NULL,
        type_ = NULL,
        response_ = NULL,
        fixpar_ = NULL,
        mats_ = NULL,
        other_data_ = NULL,
        link_ = NULL,
        invlink_ = NULL,
        coeff_fe_ = NULL,
        coeff_re_ = NULL,
        lambda_ = NULL,
        rho_ = NULL,
        terms_ = NULL,
        tmb_obj_ = NULL,
        tmb_obj_joint_= NULL,
        out_ = NULL,
        tmb_rep_ = NULL
    )
)
