
#include <TMB.hpp>
#include "nllk_sde.hpp"
#include "nllk_ctcrw.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
    //======//
    // DATA //
    //======//
    DATA_STRING(type); // Model type for each response variable
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IVECTOR(ncol_re); // Number of columns of S and X_re for each random effect
    DATA_MATRIX(a0); // Initial states for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter
    
    // Number of observations
    int n = obs.rows();
    // Time intervals (needs to be of length n for CTCRW)
    vector<Type> dtimes(n);
    for(int i = 0; i < n-1; i++)
        dtimes(i) = times(i+1) - times(i);
    dtimes(n-1) = 1;
    
    //============//
    // PARAMETERS //
    //============//
    PARAMETER_VECTOR(coeff_fe); // Fixed effect parameters
    PARAMETER_VECTOR(log_lambda); // Smoothness parameters
    PARAMETER_VECTOR(coeff_re); // Random effect parameters
    
    // Derived parameters (linear predictors)
    vector<Type> par_vec = X_fe * coeff_fe + X_re * coeff_re;
    matrix<Type> par_mat(n, par_vec.size()/n);
    for(int i = 0; i < par_mat.cols(); i++) {
        // Matrix with one row for each time step and
        // one column for each parameter
        par_mat.col(i) = par_vec.segment(i*n, n);
    }
    
    //============//
    // Likelihood //
    //============//
    Type nllk = 0;
    if (type == "BM" || type == "OU") {
        nllk = nllk_sde(ID, dtimes, obs, par_mat, type);
    } else if (type == "CTCRW") {
        nllk = nllk_ctcrw(ID, dtimes, obs, a0, P0, par_mat);
    } else {
        error ("Unknown SDE type");
    }
    
    //===================//
    // Smoothing penalty //
    // ===================//
    // Are there random effects?
    if(ncol_re(0) > 0) {
        // Index in matrix S
        int S_start = 0;
        
        // Loop over smooths
        for(int i = 0; i < ncol_re.size(); i++) {
            // Size of penalty matrix for this smooth
            int Sn = ncol_re(i);
            
            // Penalty matrix for this smooth
            Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);
            
            // Coefficients for this smooth
            vector<Type> this_coeff_re = coeff_re.segment(S_start, Sn);
            
            // Add penalty
            nllk = nllk -
                Type(0.5) * Sn * log_lambda(i) +
                Type(0.5) * exp(log_lambda(i)) * 
                density::GMRF(this_S).Quadform(this_coeff_re);
            
            // Increase index
            S_start = S_start + Sn;
        }
    }
    
    return nllk;
}
