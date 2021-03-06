
#' Create block diagonal matrix (safe version)
#' 
#' This version of bdiag checks whether the matrices passed as
#' arguments are NULL. This avoids errors that would arise if
#' using bdiag directly.
#' 
#' @param ... Matrix or list of matrices (only the first argument is used)
#' 
#' @return Block diagonal matrix
#' 
#' @importFrom Matrix bdiag
bdiag_check <- function(...) {
    # Only use first argument (matrix or list of matrices)
    args <- list(...)[[1]]
    
    # Check which matrices are non-NULL
    # (Needs two conditions to keep non-empty vectors which don't have
    # dimensions, and also keep empty matrices with set dimensions)
    check <- sapply(args, function(arg) {
        !is.null(dim(arg)) | length(arg) > 0
    })
    
    if(length(check) == 0)
        return(NULL)
    else
        return(bdiag(args[check]))
}

#' Grid of covariates
#' 
#' @param var Name of variable
#' @param data Data frame containing the covariates
#' @param covs Optional data frame with a single row and one column
#' for each covariate, giving the values that should be used. If this is
#' not specified, the mean value is used for numeric variables, and the
#' first level for factor variables.
#' @param formulas List of formulas used in the model
#' 
#' @return Data frame of covariates, with 'var' defined over a grid,
#' and other covariates fixed to their mean (numeric) or first level
#' (factor).
cov_grid <- function(var, data, covs = NULL, formulas) {
    # Get covariate names
    var_names <- unique(rapply(formulas, all.vars))
    
    # pi might appear in the formulas (e.g. used in periodic terms),
    # in which case it needs to be added to the data frame
    if(any(var_names == "pi")) {
        data$pi <- pi
    }
    
    # Get data frame of covariates
    all_vars <- data[, var_names, drop = FALSE]
    
    # Grid of covariate
    if(is.factor(all_vars[, var])) {
        n_grid <- length(unique(all_vars[, var]))
        grid <- unique(all_vars[, var])
    } else {
        n_grid <- 1e3
        grid <- seq(min(all_vars[, var]), max(all_vars[, var]), length = n_grid)
    }
    
    # New data frame for covariate grid
    new_data <- matrix(NA, nrow = n_grid, ncol = ncol(all_vars))
    colnames(new_data) <- colnames(all_vars)
    new_data <- as.data.frame(new_data)
    new_data[, var] <- grid
    
    # Select value for other covariates
    covs_list <- lapply(seq_along(all_vars), function(i) {
        if(!is.null(covs[[var_names[i]]])) {
            # Set to user-provided value if possible
            return(covs[[var_names[i]]])
        } else {
            # No user-provided value
            col <- all_vars[, i]
            if(is.numeric(col)) {
                # If numeric, use mean value
                return(mean(col, na.rm = TRUE)) 
            } else {
                # If factor, use first factor level
                return(unique(col)[1])
            }
        }
    })
    covs <- as.data.frame(covs_list)
    colnames(covs) <- colnames(all_vars)
    
    # Fill columns for other covariates
    for(var_name in colnames(new_data)) {
        if(var_name != var)
            new_data[, var_name] <- covs[, var_name]
    }
    
    return(new_data)
}

#' Compute effective degrees of freedom for a smooth term (i.e., GAM)
#'
#' This function is adapted from Dave Miller's code in the
#' package CTMCdive
#'
#' @param X_re Design matrix for random effects 
#' @param S Smoothing matrix 
#' @param lambda Smoothing parameter
#'
#' @return Trace of F = (Xt X + sp*S)^-1 Xt X
edf_smooth <- function(X_re, S, lambda) {
    # EDF = 0 if no smooth term
    if(length(lambda) == 0) {
        return(0)
    }
    
    # Duplicate lambda enough times
    lambda <- rep(lambda, nrow(S))
    
    # Calculate lambda * S
    Sbig <- S * lambda
    
    # Calculate the hat matrix
    XtX <- t(X_re) %*% X_re
    Fi <- solve(XtX + Sbig)
    F <- Fi %*% XtX
    
    # Return the trace
    return(sum(diag(F)))
}

#' logLik function for SDE objects
#' 
#' This function makes it possible to call generic R methods such
#' as AIC and BIC on SDE objects.
#' 
#' @param object SDE model object
#' @param ... For compatibility with S3 method
#' 
#' @return Maximum log-likelihood value for the model, with attributes
#' \code{df} (degrees of freedom) and \code{nobs} (number of observations)
#' 
#' @export
logLik.SDE <- function(object, ...) {
    val <- -object$out()$value 
    attributes(val)$df <- object$edf()
    attributes(val)$nobs <- nrow(object$data())
    class(val) <- "logLik"
    return(val)
}
