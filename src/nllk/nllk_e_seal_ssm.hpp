#ifndef _ESEAL_SSM_
#define _ESEAL_SSM_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Make T matrix for Kalman filter
//' 
//' @param r Drift parameter
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeT_eseal_ssm(Type r, Type dt) {
    matrix<Type> T(2, 2);
    T.setZero();
    T(0, 0) = 1;
    T(1, 0) = r * dt;
    T(1, 1) = 1;
    return T;
}

//' Make Q matrix for Kalman filter
//' 
//' @param s Diffusion parameter
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeQ_eseal_ssm(Type s, Type dt) {
    matrix<Type> Q(2, 2);
    Q.setZero();
    Q(1,1) = s * s * dt;
    return Q;
}

//' Make Z matrix for Kalman filter
//' 
//' @param a1 Parameter alpha_1 of SSM
//' @param a2 Parameter alpha_2 of SSM
//' @param m Non-lipid tissue mass
template<class Type>
matrix<Type> makeZ_eseal_ssm(Type a1, Type a2, Type m) {
    matrix<Type> Z(1, 2);
    Z(0, 0) = a1;
    Z(0, 1) = a2 / m;
    return Z;
}

//' Make H matrix for Kalman filter
//' 
//' @param tau Parameter tau of SSM
//' @param h Number of daily drift dives
template<class Type>
matrix<Type> makeH_eseal_ssm(Type tau, Type h) {
    matrix<Type> H(1, 1);
    H(0, 0) = tau * tau / h;
    return H;
}

//' Penalised negative log-likelihood for state-space model
//' of elephant seal body condition
template <class Type>
Type nllk_eseal_ssm(objective_function<Type>* obj) {
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IVECTOR(ncol_re); // Number of columns of S and X_re for each random effect
    DATA_MATRIX(a0); // Initial state estimate for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter
    DATA_VECTOR(h); // Number of daily drift dives
    DATA_VECTOR(m); // Non-lipid tissue mass
    
    // Number of observations
    int n = obs.rows();
    
    // Time intervals
    vector<Type> dtimes(n);
    for(int i = 0; i < n-1; i++)
        dtimes(i) = times(i+1) - times(i);
    dtimes(n-1) = 1;
    
    //============//
    // PARAMETERS //
    //============//
    // Other parameters of the SSM
    PARAMETER(log_tau);
    Type a1 = -0.6; // Values taken from Schick et al.
    Type a2 = 1.2; // (could not both be recovered in simulations)
    
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
    
    // Parameters of latent process
    vector<Type> r = par_mat.col(0).array();
    vector<Type> s = exp(par_mat.col(1).array());
    
    // Type a2 = exp(log_a2);
    Type tau = exp(log_tau);
    
    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> F(1, 1);
    F.setZero();
    matrix<Type> K(2, 1);
    K.setZero();
    matrix<Type> L(2, 2);
    L.setZero();
    vector<Type> u(1);
    u.setZero();
    Type detF = 0;
    
    // Initial state mean
    vector<Type> aest(2);
    aest = a0.row(0);
    int k = 1; // Counter for track (to initialise the state vector)
    // Initial state covariance matrix
    matrix<Type> Pest(2, 2);
    Pest = P0;
    
    // Kalman filter iterations
    Type llk = 0;
    for(int i = 1; i < n; i++) {
        if(ID(i) != ID(i-1)) {
            // If first location of track, re-initialise state vector
            aest = a0.row(k);
            k = k + 1;
            Pest = P0;
        } else {
            // Compute Kalman filter matrices
            matrix<Type> Z = makeZ_eseal_ssm(a1, a2, m(i));
            matrix<Type> H = makeH_eseal_ssm(tau, h(i));
            matrix<Type> T = makeT_eseal_ssm(r(i), dtimes(i));
            matrix<Type> Q = makeQ_eseal_ssm(s(i), dtimes(i));
            
            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest;
                Pest = T * Pest * T.transpose() + Q;
            } else {
                // Measurement residual
                vector<Type> obsrow =  obs.row(i).transpose();
                u = obsrow - Z * aest;
                // Residual covariance
                F = Z * Pest * Z.transpose() + H;
                detF = F(0,0);
                
                if(detF<=0) {
                    aest = T * aest;
                    Pest = T * Pest * T.transpose() + Q;
                } else {
                    // Update log-likelihood
                    matrix<Type> FinvT = F.inverse().transpose();
                    vector<Type> FinvTu = FinvT * u;
                    Type uFu = u(0) * FinvTu(0);
                    llk = llk - (log(detF) + uFu)/2;
                    // Kalman gain
                    K = T * Pest * Z.transpose() * F.inverse();
                    // Update state estimate
                    aest = T * aest + K * u;
                    // Update estimate covariance
                    L = T - K * Z;
                    Pest = T * Pest * L.transpose() + Q;
                }
            }
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
