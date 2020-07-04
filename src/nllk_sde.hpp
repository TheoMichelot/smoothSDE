#ifndef _SDE_
#define _SDE_

#include "tr_dens.hpp"

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 

//' Negative log-likelihood for SDE model
//' 
//' @param ID Vector of time series IDs
//' @param dtimes Vector of time intervals
//' @param obs Matrix of observations (two columns: x and y)
//' @param par_mat Matrix of parameters (two columns: beta and sigma)
//' @param type SDE type ("BM", "OU", ...)
//' 
//' @return Negative log-likelihood (unpenalized)
template <class Type>
Type nllk_sde(vector<Type> ID, vector<Type> dtimes, 
              matrix<Type> obs, matrix<Type> par_mat,
              std::string type) {
    // Number of observations
    int n = obs.rows();
    
    Type llk = 0;
    // Loop over observations
    for(int i = 1; i < n; i ++) {
        // No contribution if first observation of the track
        if(ID(i-1) == ID(i)) {
            llk = llk + tr_dens<Type>(obs(i, 0), obs(i-1, 0), dtimes(i-1), 
                                      par_mat.row(i-1), true, type);
        }
    }
    
    return -llk;
}

#endif
