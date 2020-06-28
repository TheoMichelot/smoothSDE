
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
        #' @description Make model matrices
        make_mats() = function() {
            return(NULL)
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
