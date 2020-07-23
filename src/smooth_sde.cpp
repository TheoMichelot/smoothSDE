
#include <TMB.hpp>
#include "nllk/nllk_sde.hpp"
#include "nllk/nllk_ctcrw.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
    // SDE type
    DATA_STRING(type);
    
    if (type == "BM" || type == "OU") {
        return nllk_sde(this);
    } else if (type == "CTCRW") {
        return nllk_ctcrw(this);
    } else {
        error ("Unknown SDE type");
    }
    return 0;
}
