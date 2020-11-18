
#ifndef _SDE_
#define _SDE_

#include "tr_dens.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Penalised negative log-likelihood for SDE model
template <class Type>
Type nllk_sde(objective_function<Type>* obj) { 
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
    DATA_VECTOR(other_data); // Optional extra data needed to evaluate the likelihood
    DATA_VECTOR(t_decay);
    DATA_IVECTOR(col_decay);
    
    // Number of observations
    int n = obs.rows();
    // Time intervals
    vector<Type> dtimes = diff(times);
    
    //============//
    // PARAMETERS //
    //============//
    PARAMETER_VECTOR(coeff_fe); // Fixed effect parameters
    PARAMETER_VECTOR(log_lambda); // Smoothness parameters
    PARAMETER(log_decay);
    PARAMETER_VECTOR(coeff_re); // Random effect parameters
    
    Type decay_rate = exp(log_decay);
    
    matrix<Type> X_re_copy = X_re;
    
    if(t_decay.size() > 1) {
        for(int i = 0; i < col_decay.size(); i++) {
            int i_col = col_decay(i) - 1;
            vector<Type> decay = exp(-decay_rate * t_decay);
            vector<Type> X_col = X_re.col(i_col);
            vector<Type> X_decay = X_col * decay;
            X_re_copy.col(i_col) = X_decay;
        }
    }
    
    // Derived parameters (linear predictors)
    vector<Type> par_vec = X_fe * coeff_fe + X_re_copy * coeff_re;
    matrix<Type> par_mat(n, par_vec.size()/n);
    for(int i = 0; i < par_mat.cols(); i++) {
        // Matrix with one row for each time step and
        // one column for each parameter
        par_mat.col(i) = par_vec.segment(i*n, n);
    }
    
    //============//
    // Likelihood //
    //============//
    // Initialise log-likelihood
    Type llk = 0;
    // Loop over observations
    for(int i = 1; i < n; i ++) {
        // No contribution if first observation of the track
        if(ID(i-1) == ID(i)) {
            llk = llk + tr_dens<Type>(obs(i, 0), obs(i-1, 0), dtimes(i-1), 
                                      par_mat.row(i-1), true, type, other_data);
        }
    }
    
    //===================//
    // Smoothing penalty //
    // ===================//
    Type nllk = -llk;
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

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif