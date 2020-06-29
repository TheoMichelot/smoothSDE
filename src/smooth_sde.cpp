
#include <TMB.hpp>
#include "TMB/std_sde.hpp"
#include "TMB/ctcrw.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_STRING(type);
    if (type == "BM" || type == "OU") {
        return std_sde(this);
    } else if (type == "CTCRW") {
        return ctcrw(this);
    } else {
        error ("Unknown SDE type");
    }
    return 0;
}
