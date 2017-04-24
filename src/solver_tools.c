
#include"solver_tools.h"

double L2norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity){

    double L2norm_=0.0;

    register int i; int j;

    for(i=0; i<Nelem_; i++)
        for(j=0; j<Ndof_per_elem; j++)
            L2norm_ += pow(quantity[i][j],2);

    L2norm_ = sqrt(L2norm_/Ndof_per_elem);

    return L2norm_;
}

double L2norm_perdof(const int& Nelem_, double* quantity){

    double L2norm_=0.0;

    register int i;

    for(i=0; i<Nelem_; i++)
            L2norm_ += pow(quantity[i],2);

    L2norm_ = sqrt(L2norm_/Nelem_);

    return L2norm_;
}


