
#include <TMB.hpp>
#include "TMB/std_sde.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_IVECTOR(type);
    if (type(0) == 1 || type(0) == 2) {
        return std_sde(this);
    } else {
        error ("Unknown SDE type");
    }
    return 0;
}
