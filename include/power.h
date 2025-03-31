#ifndef POWER_H
#define POWER_H

#include "matrix.h"

void print_XAX(
    const double *A, const double *X, int n, int k, int nb, char *msg
);

void compute_eigv_shift(
    SymBandMatrix *K,
    SymBandMatrix *M,
    double *eigw,
    double *eigv,
    double rtol,
    int max_iter
);

void compute_eigvs_deflation(
    SymBandMatrix *K,
    SymBandMatrix *M,
    double *eigws,
    double *eigvs,
    int nb,
    double rtol,
    int max_iter
);

#endif