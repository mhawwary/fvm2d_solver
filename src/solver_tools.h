#pragma once

#include"general_tools.h"

double L1norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol);

double L1norm_perdof(const int& Nelem_, double* quantity);

double L2norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol);

double L2norm_perdof(const int& Nelem_, double* quantity);
