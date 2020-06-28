
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
        
        #' @description Setup model
        setup = function() {
            return(NULL)
        },
        
        #' @description Fit model 
        fit = function() {
            return(NULL)
        }
    ),
    
    private = list(
        formulas_ = NULL,
        data_ = NULL,
        type_ = NULL
    )
)
