
#ifndef _TR_DENS_
#define _TR_DENS_

template<class Type>
Type tr_dens(Type Z1, Type Z0, Type dtimes, vector<Type> par, 
             bool give_log, std::string type, vector<Type> other_data) {
    Type res = 0;
    Type mean = 0;
    Type sd = 1;
    
    if(type == "BM") {
        // Brownian motion: dZ_t = mu(t) dt + sigma(t) dW_t
        // where par = (mu, log(sigma))
        mean = Z0 + par(0) * dtimes;
        sd = exp(par(1)) * sqrt(dtimes);
        res = dnorm(Z1, mean, sd, give_log);
    } else if(type == "BM-t") {
        // Brownian motion with t-distributed noise
        Type df = other_data(0); // number of degrees of freedom
        mean = par(0) * dtimes;
        sd = exp(par(1)) * sqrt(dtimes);
        Type scale = sd/sqrt(df/(df-2));
        res = dt((Z1-Z0-mean)/scale, df, give_log) - log(scale);
    } else if(type == "OU") {
        // Ornstein-Uhlenbeck: dZ_t = beta(t) (mu(t) - Z_t) dt + sigma(t) dW_t
        // where par = mu, log(beta), log(sigma)
        mean = par(0) + exp(- exp(par(1)) * dtimes) * (Z0 - par(0));
        sd = exp(par(2))/sqrt(2 * exp(par(1))) * 
            sqrt(1 - exp(-2 * exp(par(1)) * dtimes));
        res = dnorm(Z1, mean, sd, give_log);
    }
    
    return(res);
}

#endif