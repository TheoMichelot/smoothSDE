
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
        },
        
        ###############
        ## Accessors ##
        ###############
        #' @description Return formulas of SDE object
        formulas = function() {return(private$formulas_)},
        
        #' @description Return data of SDE object
        data = function() {return(private$data_)},
        
        #' @description Return type of SDE object
        type = function() {return(private$type_)},
        
        #' @description Return type of SDE object as integer code
        type_code = function() {
            switch (self$type(),
                    "BM" = 1,
                    "OU" = 2)
        },
        
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
        make_mats = function(new_data = NULL) {
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
        
        #' @description TMB setup
        #'  
        #' This creates an attribute \code{tmb_obj}, which can be used to 
        #' evaluate the negative log-likelihood function.
        #' 
        #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
        setup = function(silent = TRUE) {
            # Create model matrices
            mats <- self$make_mats()
            X_fe <- mats$X_fe
            X_re <- mats$X_re
            S <- mats$S
            ncol_re <- mats$ncol_re
            
            # Format initial parameters for TMB
            tmb_par <- list(coeff_fe = rep(0, ncol(X_fe)),
                            coeff_re = 0,
                            log_lambda = 0)
            
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
            private$tmb_rep_ <- sdreport(private$tmb_obj_)
            par_list <- as.list(private$tmb_rep_, "Estimate")
        }
    ),
    
    private = list(
        formulas_ = NULL,
        data_ = NULL,
        type_ = NULL,
        tmb_obj_ = NULL,
        fit_ = NULL,
        tmb_rep_ = NULL
    )
)
