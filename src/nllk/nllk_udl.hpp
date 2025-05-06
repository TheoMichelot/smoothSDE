#ifndef _UDL_
#define _UDL_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Make T matrix for Kalman filter
 //' @param gamma Parameter gamma of UDL process
 //' @param dt Length of time interval
 //' @param n_dim Number of dimensions of UDL process
 template<class Type>
 matrix<Type> makeT_udl(Type gamma, Type dt, int n_dim) {
     matrix<Type> T(2*n_dim, 2*n_dim);
     T.setZero();
     for(int i = 0; i < n_dim; i++) {
         T(2*i, 2*i) = 1;
         T(2*i, 2*i + 1) = (1 - exp(- gamma * dt)) / gamma;
         T(2*i + 1, 2*i + 1) = exp(- gamma * dt);
     }
     return T;
 }
 
 //' Make B matrix for Kalman filter
 //' 
 //' @param gamma Parameter gamma of UDL process
 //' @param sigma Parameter sigma of UDL process
 //' @param dt Length of time interval
 //' @param n_dim Number of dimensions of UDL process
 template<class Type>
 matrix<Type> makeB_udl(Type gamma, Type sigma, Type dt, int n_dim) {
     matrix<Type> B(2*n_dim, n_dim);
     B.setZero();
     for(int i = 0; i < n_dim; i++) {
         B(2*i, i) = sigma*sigma / gamma * (dt - 1/gamma * (1 - exp(-gamma * dt)));
         B(2*i + 1, i) = sigma*sigma / gamma * (1 - exp(-gamma * dt));        
     }
     return B;
 }
 
 //' Make Q matrix for Kalman filter
 //' 
 //' @param gamma Parameter gamma of UDL process
 //' @param sigma Parameter sigma of UDL process
 //' @param dt Length of time interval
 //' @param n_dim Number of dimensions of UDL process
 template<class Type>
 matrix<Type> makeQ_udl(Type gamma, Type sigma, Type dt, int n_dim) {
     matrix<Type> Q(2*n_dim, 2*n_dim);
     Q.setZero();
     Type sig2 = sigma * sigma;
     Type gamma2 = gamma * gamma;
     for(int i = 0; i < n_dim; i++) {
         Q(2*i, 2*i) = sig2 * (2 * dt / gamma - exp(-2 * gamma * dt)/gamma2 -
             3 / gamma2 + 4 * exp(- gamma * dt) / gamma2);
         Q(2*i, 2*i + 1) = sig2 / gamma * 
             (1 - 2 * exp(- gamma * dt) + exp(- 2 * gamma * dt));
         Q(2*i + 1, 2*i) = Q(2*i, 2*i + 1);
         Q(2*i + 1, 2*i + 1) = sig2 * (1 - exp(-2 * gamma * dt));
     }
     return Q;
 }
 
 //' Negative log-likelihood for UDL process
 template<class Type>
 Type nllk_udl(objective_function<Type>* obj) {
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
     DATA_ARRAY(cov_grad);
     
     DATA_ARRAY(H_array); // Covariance matrices for observation error
     
     // Prior information 
     DATA_MATRIX(coeff_fe_prior); // means, sds for prior on fixed effects
     DATA_MATRIX(log_lambda_prior); // means, sds for prior on smoothing parameters
     
     // Number of observations
     int n = obs.rows();
     
     // Number of dimensions
     int n_dim = obs.cols();
     
     // Time intervals (needs to be of length n)
     vector<Type> dtimes(n);
     for(int i = 0; i < n-1; i++) {
         if(ID(i+1) == ID(i)) {
             dtimes(i) = times(i+1) - times(i);
         } else {
             // Use last time interval twice
             dtimes(i) = times(i) - times(i-1);
         }
     }
     dtimes(n-1) = dtimes(n-2);
     
     //============//
     // PARAMETERS //
     //============//
     // // SD of measurement error
     // PARAMETER(log_sigma_obs);
     // Type sigma_obs = exp(log_sigma_obs);
     
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
     
     // Parameters of UDL process
     int n_cov = cov_grad.cols();
     vector<Type> gamma = exp(par_mat.col(0).array());
     vector<Type> sigma = exp(par_mat.col(1).array());
     matrix<Type> beta = par_mat.block(0, 2, n, n_cov + 2).array();
     
     // Gradient of stationary distribution
     matrix<Type> h(n, n_dim);
     h.setZero();
     for(int i = 0; i < n; i++) {
         for(int j = 0; j < n_dim; j++) {
             for(int k = 0; k < n_cov; k++) {
                 h(i, j) = h(i, j) + beta(i, k) * cov_grad(i, j, k);                 
             }
         }
     }
     
     //================================//
     // Likelihood using Kalman filter //
     //================================//
     // Define all matrices and vectors needed below
     matrix<Type> Z(n_dim, 2*n_dim);
     Z.setZero();
     for(int i = 0; i < n_dim; i++) {
         Z(i, 2*i) = 1;
     }
     // matrix<Type> H = makeH_udl(sigma_obs, n_dim);
     matrix<Type> T(2*n_dim, 2*n_dim);
     matrix<Type> Q(2*n_dim, 2*n_dim);
     matrix<Type> B(2*n_dim, n_dim);
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
     matrix<Type> aest_all(n, 2*n_dim);
     aest_all.setZero();
     aest_all.row(0) = aest;
     matrix<Type> residuals(n, n_dim);
     residuals.setZero();
     for(int i = 1; i < n; i++) {
         if(ID(i) != ID(i-1)) {
             // If first location of track, re-initialise state vector
             aest = a0.row(k);
             k = k + 1;
             Pest = P0;
         } else {
             // Compute Kalman filter matrices
             // if(H_array.size() > 1) {
             //     H = H_array.col(i).matrix();
             // }
             T = makeT_udl(gamma(i), dtimes(i), n_dim);
             Q = makeQ_udl(gamma(i), sigma(i), dtimes(i), n_dim);
             B = makeB_udl(gamma(i), sigma(i), dtimes(i), n_dim);
             
             vector<Type> h_i = h.row(i-1).transpose();
             vector<Type> B_times_h = B * h_i;
             
             if(R_IsNA(asDouble(obs(i,0)))) {
                 // If missing observation
                 aest = T * aest + B_times_h;
                 Pest = T * Pest * T.transpose() + Q;
             } else {
                 // Measurement residual
                 vector<Type> obsrow =  obs.row(i).transpose();
                 u = obsrow - Z * aest;
                 residuals.row(i) = u;
                 // Residual covariance
                 F = Z * Pest * Z.transpose(); //+ H;
                 detF = det(F);
                 
                 if(detF <= 0) {
                     aest = T * aest + B_times_h;
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
                     aest = T * aest + K * u + B_times_h;
                     // Update estimate covariance
                     L = T - K * Z;
                     Pest = T * Pest * L.transpose() + Q;
                 }
             }
         }        
         
         aest_all.row(i) = aest;
     }
     
     REPORT(residuals)
     
     //========//
     // Priors //
     //========//
     // Fixed effects
     for (int i = 0; i < coeff_fe.size(); i++) {
         if (!R_IsNA(asDouble(coeff_fe_prior(i, 0)))) {
             llk = llk + dnorm(coeff_fe(i), coeff_fe_prior(i, 0), coeff_fe_prior(i, 1), 1.0); 
         }
     }
     // Smoothing parameters
     if (ncol_re(0) > 0) {
         for (int i = 0; i < log_lambda.size(); i++) {
             if (!R_IsNA(asDouble(log_lambda_prior(i, 0)))) {
                 llk = llk + dnorm(log_lambda(i), log_lambda_prior(i, 0), log_lambda_prior(i, 1), 1.0); 
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