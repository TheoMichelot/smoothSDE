
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
        #' @param formulas Nested list of formulas for model
        #' parameters
        #' @param data Data frame with covariates, response variables,
        #' time, and ID
        #' @param models List of model names for observed processes
        #' 
        #' @return A new Dist object
        initialize = function(formulas, data, models) {
            private$formulas_ <- formulas
            private$data_ <- data
            private$models_ <- models
        },
        
        ###############
        ## Accessors ##
        ###############
        #' @description Return name of Dist object
        formulas = function() {return(private$formulas_)},
        
        #' @description Return pdf of Dist object
        data = function() {return(private$data_)},
        
        #' @description Return random generator function of Dist object
        models = function() {return(private$models_)},
        
        ###################
        ## Other methods ##
        ###################
        #' @description Fit model 
        fit = function() {
            return(NULL)
        }
    ),
    
    private = list(
        formulas_ = NULL,
        data_ = NULL,
        models_ = NULL
    )
)
