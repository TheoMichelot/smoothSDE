#ifndef _CTCRW_
#define _CTCRW_

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

//' Negative log-likelihood for CTCRW
//' 
//' This function was inspired by the source code of the package
//' crawl, authored by Devin Johnson and Josh London
//' 
//' @param ID Vector of time series IDs
//' @param dtimes Vector of time intervals
//' @param obs Matrix of observations (two columns: x and y)
//' @param a0 Matrix of initial state vectors for Kalman filter
//' @param P0 Initial state covariance matrix for Kalman filter
//' @param par_mat Matrix of parameters (two columns: beta and sigma)
//' 
//' @return Negative log-likelihood (unpenalized)
template <class Type>
Type nllk_ctcrw(vector<Type> ID, vector<Type> dtimes, 
                matrix<Type> obs, matrix<Type> a0,
                matrix<Type> P0, matrix<Type> par_mat) {
    // Number of observations
    int n = obs.rows();
    
    // Parameters of velocity process
    vector<Type> beta = exp(par_mat.col(0).array());
    vector<Type> sigma = exp(par_mat.col(1).array());
    
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
    aest = a0.row(1);
    // Initial state covariance matrix
    matrix<Type> Pest(4,4);
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
    
    return -llk;
}

#endif
