
#include"solver_tools.h"

double L1norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol){

    double L1norm_=0.0;

    register int i; int j;

    for(i=0; i<Nelem_; i++)
        for(j=0; j<Ndof_per_elem; j++)
            L1norm_ += fabs(quantity[i][j]*Vol[i]);

    L1norm_ = L1norm_/(Ndof_per_elem*Nelem_);

    return L1norm_;
}

double L1norm_perdof(const int& Nelem_, double* quantity){

    double L1norm_=0.0;

    register int i;

    for(i=0; i<Nelem_; i++)
            L1norm_ += fabs(quantity[i]);

    L1norm_ = L1norm_/Nelem_;

    return L1norm_;
}

double L2norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol){

    double L2norm_=0.0;

    register int i; int j;

    for(i=0; i<Nelem_; i++)
        for(j=0; j<Ndof_per_elem; j++)
            L2norm_ += pow(quantity[i][j]*Vol[i],2);

    L2norm_ = sqrt(L2norm_/(Ndof_per_elem*Nelem_));

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


