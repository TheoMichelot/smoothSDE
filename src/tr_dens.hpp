
#ifndef _TR_DENS_
#define _TR_DENS_

template<class Type>
Type tr_dens(Type Z1, Type Z0, Type dtimes, vector<Type> par, 
             bool log, std::string type) {
    Type res = 0;
    Type mean = 0;
    Type sd = 1;
    
    if(type == "BM") {
        // Brownian motion: dZ_t = mu(t) dt + sigma(t) dW_t
        // where par = (mu, log(sigma))
        mean = Z0 + par(0) * dtimes;
        sd = exp(par(1)) * sqrt(dtimes);
        res = dnorm(Z1, mean, sd, log);
    } else if(type == "OU") {
        // Ornstein-Uhlenbeck: dZ_t = beta(t) (mu(t) - Z_t) dt + sigma(t) dW_t
        // where par = mu, log(beta), log(sigma)
        mean = par(0) + exp(- exp(par(1)) * dtimes) * (Z0 - par(0));
        sd = exp(par(2))/(2 * exp(par(1))) * (1 - exp(-2 * exp(par(1)) * dtimes));
        res = dnorm(Z1, mean, sd, log);
    }
    
    return(res);
}

#endif