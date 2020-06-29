
#ifndef _CTCRW_
#define _CTCRW_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Make T matrix for Kalman filter
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeT(Type beta, Type dt) {
    matrix<Type> T(4,4);
    T.setZero();
    
    T(0,0) = 1;
    T(2,2) = 1;
    T(0,1) = (1-exp(-beta*dt))/beta;
    T(1,1) = exp(-beta*dt);
    T(2,3) = (1-exp(-beta*dt))/beta;
    T(3,3) = exp(-beta*dt);
    
    return T;
}

//' Make Q matrix for Kalman filter
//' 
//' @param beta Parameter beta of OU velocity process
//' @param sigma Parameter sigma of OU velocity process
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeQ(Type beta, Type sigma, Type dt) {
    matrix<Type> Q(4,4);
    Q.setZero();
    
    Q(0,0) = (sigma/beta)*(sigma/beta)*(dt - 2/beta*(1-exp(-beta*dt)) + 
        1/(2*beta)*(1-exp(-2*beta*dt)));
    Q(1,1) = sigma*sigma/(2*beta) * (1-exp(-2*beta*dt));
    Q(0,1) = sigma*sigma/(2*beta*beta) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt));
    Q(1,0) = Q(0,1);
    
    Q(2,2) = (sigma/beta)*(sigma/beta)*(dt - 2/beta*(1-exp(-beta*dt)) + 
        1/(2*beta)*(1-exp(-2*beta*dt)));
    Q(3,3) = sigma*sigma/(2*beta) * (1-exp(-2*beta*dt));
    Q(2,3) = sigma*sigma/(2*beta*beta) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt));
    Q(3,2) = Q(2,3);
    
    return Q;
}

template <class Type>
Type ctcrw(objective_function<Type>* obj) {
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
    DATA_VECTOR(a0); // Initial state estimate for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter
    
    // Number of observations
    int n = obs.rows();
    // Time intervals (needs to be of length n)
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
    
    // Parameters of velocity process
    vector<Type> beta = exp(par_mat.col(0).array());
    vector<Type> sigma = exp(par_mat.col(1).array());

    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> Z(2, 4);
    Z.setZero();
    Z(0,0) = 1;
    Z(1,2) = 1;
    matrix<Type> H(2,2);
    H.setZero(); // for now, no measurement error
    matrix<Type> T(4,4);
    matrix<Type> Q(4,4);
    matrix<Type> F(2,2);
    F.setZero();
    matrix<Type> K(4,2);
    K.setZero();
    matrix<Type> L(4,4);
    L.setZero();
    vector<Type> u(2);
    u.setZero();
    Type detF;
    
    // Initial state mean
    vector<Type> aest(4);
    aest = a0;
    // Initial state covariance matrix
    matrix<Type> Pest(4,4);
    Pest = P0;
    
    // Kalman filter iterations
    Type llk = 0;
    for(int i = 1; i < n; i++) {
        if(ID(i) != ID(i-1)) {
            // If first location of track, re-initialise state vector
            aest = a0;
            Pest = P0;
        } else {
            // Compute Kalman filter matrices
            matrix<Type> T = makeT(beta(i), dtimes(i));
            matrix<Type> Q = makeQ(beta(i), sigma(i), dtimes(i));
            
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
                detF = F(0,0) * F(1,1) - F(1,0) * F(0,1);
                
                if(detF<=0) {
                    aest = T * aest;
                    Pest = T * Pest * T.transpose() + Q;
                } else {
                    // Update log-likelihood
                    matrix<Type> FinvT = F.inverse().transpose();
                    vector<Type> FinvTu = FinvT * u;
                    Type uFu = u(0) * FinvTu(0) + u(1) * FinvTu(1);
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
