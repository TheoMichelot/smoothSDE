
#ifndef _TR_DENS_
#define _TR_DENS_

//' Transition density for a few basic diffusion processes
//' 
//' @param Z1 Vector of observations at time t_{i+1}
//' @param Z0 Vector of observations at time t_i
//' @param dtimes Time interval = t_{i+1} - t_i
//' @param par Vector of SDE parameters on working scale (transformed to
//' natural scale within this function)
//' @param give_log Logical. Should the log-density be returned?
//' @param type String for model type ("BM", "BM-t" or "OU")
//' @param other_data Vector for any additional data needed; e.g., the number
//' of degrees of freedom for BM-t
//' 
//' @return Transition density from Z0 to Z1
template<class Type>
Type tr_dens(vector<Type> Z1, vector<Type> Z0, Type dtimes, vector<Type> par, 
             bool give_log, std::string type, vector<Type> other_data) {
    Type res = 0;
    Type mean = 0;
    Type sd = 1;
    
    // Number of dimensions
    int n_dim = Z1.size();
    
    // Loop over dimensions to add contribution of each
    for(int i = 0; i < n_dim; i++) {
        // No contribution if either start or end observation is missing
        if(!R_IsNA(asDouble(Z0(i))) && !R_IsNA(asDouble(Z1(i)))) {
            if(type == "BM") {
                // Brownian motion: dZ_t = mu(t) dt + sigma(t) dW_t
                // where par = (mu, log(sigma))
                mean = Z0(i) + par(0) * dtimes;
                sd = exp(par(1)) * sqrt(dtimes);
                res = res + dnorm(Z1(i), mean, sd, true);
            } else if(type == "BM-t") {
                // Brownian motion with t-distributed noise
                Type df = other_data(0); // number of degrees of freedom
                mean = par(0) * dtimes;
                sd = exp(par(1)) * sqrt(dtimes);
                Type scale = sd/sqrt(df/(df-2));
                res = res + dt((Z1(i)-Z0(i)-mean)/scale, df, true) - log(scale);
            } else if(type == "OU") {
                // Ornstein-Uhlenbeck: 
                // dZ_t = 1/tau(t) (mu(t) - Z_t) dt + sqrt(2*kappa(t)/tau(t)) dW_t
                // where par = mu_1, mu_2, ..., mu_d, log(tau), log(kappa)
                mean = par(i) + exp(- dtimes/exp(par(n_dim))) * (Z0(i) - par(i));
                sd = sqrt(exp(par(n_dim+1)) * 
                    (1 - exp(-2 * dtimes / exp(par(n_dim)))));
                res = res + dnorm(Z1(i), mean, sd, true);
            }            
        }
    }
    
    if(!give_log)
        res = exp(res);
    
    return(res);
}

#endif