#ifndef _CTCRW_
#define _CTCRW_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Matrix determinant
template<class Type>
Type det(matrix<Type> M) {
    int n_dim = M.cols();
    Type det = 0;
    if(n_dim == 1) {
        det = M(0, 0);
    } else if(n_dim == 2) {
        det = M(0,0) * M(1,1) - M(1,0) * M(0,1);        
    } else {
        det = exp(atomic::logdet(M));
    }
    return det;
}

//' Make T matrix for Kalman filter
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeT_ctcrw(Type beta, Type dt, int n_dim) {
    matrix<Type> T(2*n_dim, 2*n_dim);
    T.setZero();
    for(int i = 0; i < n_dim; i++) {
        T(2*i, 2*i) = 1;
        T(2*i, 2*i + 1) = (1-exp(-beta*dt))/beta;
        T(2*i + 1, 2*i + 1) = exp(-beta*dt);
    }
    return T;
}

//' Make Q matrix for Kalman filter
//' 
//' @param beta Parameter beta of OU velocity process
//' @param sigma Parameter sigma of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeQ_ctcrw(Type beta, Type sigma, Type dt, int n_dim) {
    matrix<Type> Q(2*n_dim, 2*n_dim);
    Q.setZero();
    for(int i = 0; i < n_dim; i++) {
        Q(2*i, 2*i) = (sigma/beta)*(sigma/beta)*(dt - 2/beta*(1-exp(-beta*dt)) + 
            1/(2*beta)*(1-exp(-2*beta*dt)));
        Q(2*i, 2*i + 1) = sigma*sigma/(2*beta*beta) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt));
        Q(2*i + 1, 2*i) = Q(2*i, 2*i + 1);
        Q(2*i + 1, 2*i + 1) = sigma*sigma/(2*beta) * (1-exp(-2*beta*dt));
    }
    return Q;
}

//' Make B matrix for Kalman filter
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeB_ctcrw(Type beta, Type dt, int n_dim) {
    matrix<Type> B(2*n_dim, n_dim);
    B.setZero();
    for(int i = 0; i < n_dim; i++) {
        B(2*i, i) = 1 - exp(-beta*dt);
        B(2*i + 1, i) = dt - (1 - exp(-beta*dt))/beta;
    }
    return B;
}

//' Penalised negative log-likelihood for CTCRW
//' 
//' This function was inspired by the source code of the package
//' crawl, authored by Devin Johnson and Josh London
//' 
//' All derivations, including for the matrices T and Q defined above, are detailed 
//' in Section 6.2.2 of Michelot (2019), Stochastic models of animal movement and 
//' habitat selection. PhD thesis, University of Sheffield.
//' (etheses.whiterose.ac.uk/23688/1/TheoMichelot_PhD_thesis_April2019.pdf)
template <class Type>
Type nllk_ctcrw(objective_function<Type>* obj) {
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
    
    // Number of observations
    int n = obs.rows();
    
    // Number of dimensions
    int n_dim = obs.cols();
    
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
    matrix<Type> mu = par_mat.block(0, 0, n, n_dim).array();
    vector<Type> beta = exp(par_mat.col(n_dim).array());
    vector<Type> sigma = exp(par_mat.col(n_dim + 1).array());
    
    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> Z(n_dim, 2*n_dim);
    Z.setZero();
    for(int i = 0; i < n_dim; i++) {
        Z(i, 2*i) = 1;
    }
    matrix<Type> H(n_dim, n_dim);
    H.setZero(); // for now, no measurement error
    matrix<Type> T(2*n_dim, 2*n_dim);
    matrix<Type> Q(2*n_dim, 2*n_dim);
    matrix<Type> F(n_dim, n_dim);
    F.setZero();
    matrix<Type> K(2*n_dim, n_dim);
    K.setZero();
    matrix<Type> L(2*n_dim, 2*n_dim);
    L.setZero();
    vector<Type> u(n_dim);
    u.setZero();
    Type detF;
    
    // Initial state mean
    vector<Type> aest(2*n_dim);
    aest = a0.row(0);
    // Initial state covariance matrix
    matrix<Type> Pest(2*n_dim, 2*n_dim);
    Pest = P0;
    
    // Counter for ID (to initialise a0)
    int k = 1;
    
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
            matrix<Type> T = makeT_ctcrw(beta(i), dtimes(i), n_dim);
            matrix<Type> Q = makeQ_ctcrw(beta(i), sigma(i), dtimes(i), n_dim);
            matrix<Type> B = makeB_ctcrw(beta(i), dtimes(i), n_dim);
            
            // Mean velocity component of state update
            vector<Type> mu_i = mu.row(i).transpose();
            vector<Type> B_times_mu = B * mu_i;

            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest + B_times_mu;
                Pest = T * Pest * T.transpose() + Q;
            } else {
                // Measurement residual
                vector<Type> obsrow =  obs.row(i).transpose();
                u = obsrow - Z * aest;
                // Residual covariance
                F = Z * Pest * Z.transpose() + H;
                detF = det(F);

                if(detF <= 0) {
                    aest = T * aest;
                    Pest = T * Pest * T.transpose() + Q;
                } else {
                    // Update log-likelihood
                    matrix<Type> FinvT = F.inverse().transpose();
                    vector<Type> FinvTu = FinvT * u;
                    Type uFu = (u * FinvTu).sum();
                    llk = llk - (log(detF) + uFu)/2;
                    // Kalman gain
                    K = T * Pest * Z.transpose() * F.inverse();
                    // Update state estimate
                    aest = T * aest + K * u + B_times_mu;
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
