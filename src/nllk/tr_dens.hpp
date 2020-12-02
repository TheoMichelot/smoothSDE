
#ifndef _TR_DENS_
#define _TR_DENS_

//' Transition density for a few basic diffusion processes
//' 
//' Z1 and Z0 are vectors to allow for multivariate observations
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
            // Ornstein-Uhlenbeck: dZ_t = beta(t) (mu(t) - Z_t) dt + sigma(t) dW_t
            // where par = mu_1, mu_2, ..., mu_d, log(beta), log(sigma)
            mean = par(i) + exp(- exp(par(n_dim)) * dtimes) * (Z0(i) - par(i));
            sd = exp(par(n_dim+1))/sqrt(2 * exp(par(n_dim))) * 
                sqrt(1 - exp(-2 * exp(par(n_dim)) * dtimes));
            res = res + dnorm(Z1(i), mean, sd, true);
        }
    }
    
    if(!give_log)
        res = exp(res);
    
    return(res);
}

#endif