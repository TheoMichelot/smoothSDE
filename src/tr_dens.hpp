
#ifndef _TR_DENS_
#define _TR_DENS_

#ifndef _TMB_
#define _TMB_
#include <TMB.hpp>
#endif 

template<class Type>
Type tr_dens(Type Z1, Type Z0, Type dtimes, vector<Type> par, bool log, int type) {
    Type res = 0;
    Type mean = 0;
    Type sd = 1;
    
    if(type == 1) {
        // Brownian motion: dZ_t = mu(t) dt + sigma(t) dW_t
        // where par = (mu, sigma)
        mean = Z0 + par(0) * dtimes;
        sd = par(1) * sqrt(dtimes);
        res = dnorm(Z1, mean, sd, log);
    } else if(type == 2) {
        // Ornstein-Uhlenbeck: dZ_t = beta(t) (mu(t) - Z_t) dt + sigma(t) dW_t
        // where par = mu, beta, sigma
        mean = par(0) + exp(- par(1) * dtimes) * (Z0 - par(0));
        sd = par(2)/(2 * par(1)) * (1 - exp(-2 * par(1) * dtimes));
        res = dnorm(Z1, mean, sd, log);
    }
    
    return(res);
}

#endif