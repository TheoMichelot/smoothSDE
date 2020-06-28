
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{   
    using namespace R_inla; 
    using namespace density; 
    using namespace Eigen; 
    
    DATA_MATRIX(Z);
    PARAMETER(beta);
    
    Type nllk = 0;
    
    return nllk;
}
