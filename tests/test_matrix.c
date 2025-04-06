#define _POSIX_C_SOURCE 199309L
#include "matrix.h"
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define RWMJ CblasRowMajor
#define LOWR CblasLower
#define UNIT CblasUnit
#define NNUN CblasNonUnit
#define TRNS CblasTrans
#define NTRS CblasNoTrans

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define TIMESP(t1, t2)                                                         \
    ((double)((t2).tv_sec - (t1).tv_sec) +                                     \
     1e-9 * ((double)((t2).tv_nsec - (t1).tv_nsec)))

void random_mat(idx_t n, idx_t k, double ratio, CSRMatrix **csr_ptr) {
    k = MIN(k, n - 1);
    SymBandMatrix *K = allocate_band_sym(n, k);
    idx_t bound;
    for (idx_t i = 0; i < n; i++) {
        bound = (i < k) ? 0 : i - k;
        for (idx_t j = bound; j < i; j++) {
            idx_t col = rand() % i;
            K->a[i][j] = -(double)rand() / (double)RAND_MAX;
            if (fabs(K->a[i][j]) < 1. - ratio / k) {
                K->a[i][j] = 0.0;
            }
        }
        K->a[i][i] = (k+1) / 2. + (k+1) / 2. * (double)rand() / (double)RAND_MAX;
    }

    CSRMatrix *mat = band_to_csr_sym(K);

    write_band_sym(K, NULL, "K.txt");
    write_csr(mat, NULL, "K_csr.txt");

    free_band_sym(K);
    *csr_ptr = mat;
}

double *csr_to_full(CSRMatrix *csr, double *L) {
    for (idx_t i = 0; i < csr->n; i++) {
        for (idx_t j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++) {
            idx_t idx = i * csr->n + csr->col_idx[j];
            L[idx] = csr->data[j];
        }
    }
}

void random_vec(double *v, idx_t n) {
    for (idx_t i = 0; i < n; i++) {
        v[i] = (double)rand() / (double)RAND_MAX;
        // v[i] = round(v[i] * 100.0) / 10.;
    }
}

int test_csr_solve(idx_t n) {

    CSRMatrix *csr;
    double *L_full, err;
    random_mat(n, n - 1, 4., &csr);
    print_csr(csr);
    printf("\n");

    L_full = (double *)calloc(n * n, sizeof(double));
    csr_to_full(csr, L_full);

    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    random_vec(b, n);
    // print_vector(n, b);
    double norm_b = cblas_dnrm2(n, b, 1);

    memcpy(x, b, n * sizeof(double));
    solve_csr(csr, x, CsrNoTrans, CsrNonUnit);
    // print_vector(n, x);
    cblas_dtrmv(RWMJ, LOWR, NTRS, NNUN, n, L_full, n, x, 1);
    cblas_daxpy(n, -1., b, 1, x, 1);
    err = cblas_dnrm2(n, x, 1) / norm_b;
    printf("Solve  Lower x = b | err = %.3le\n", err);
    
    memcpy(x, b, n * sizeof(double));
    solve_csr(csr, x, CsrNoTrans, CsrUnit);
    cblas_dtrmv(RWMJ, LOWR, NTRS, UNIT, n, L_full, n, x, 1);
    cblas_daxpy(n, -1., b, 1, x, 1);
    err = cblas_dnrm2(n, x, 1) / norm_b;
    printf("Solve  Lunit x = b | err = %.3le\n", err);

    memcpy(x, b, n * sizeof(double));
    solve_csr(csr, x, CsrTrans, CsrNonUnit);
    cblas_dtrmv(RWMJ, LOWR, TRNS, NNUN, n, L_full, n, x, 1);
    cblas_daxpy(n, -1., b, 1, x, 1);
    err = cblas_dnrm2(n, x, 1) / norm_b;
    printf("Solve  Upper x = b | err = %.3le\n", err);

    memcpy(x, b, n * sizeof(double));
    solve_csr(csr, x, CsrTrans, CsrUnit);
    cblas_dtrmv(RWMJ, LOWR, TRNS, UNIT, n, L_full, n, x, 1);
    cblas_daxpy(n, -1., b, 1, x, 1);
    err = cblas_dnrm2(n, x, 1) / norm_b;
    printf("Solve  Uunit x = b | err = %.3le\n", err);

    free(L_full);
    free(b);
    free(x);
    free_csr(csr);
    return 0;
}

void test_ilu(idx_t n, idx_t k, double ratio) {
    CSRMatrix *mat;
    random_mat(n, k, ratio, &mat);
    printf("Matrix A\n");
    print_csr(mat);

    CSRMatrix *ilu0 = ilu0_symmetric(mat);
    printf("Matrix ILU0\n");
    print_csr(ilu0);
    
    CSRMatrix *ilu1 = ilu1_symmetric(mat);
    printf("Matrix ILU1\n");
    print_csr(ilu1);

    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    double *y = (double *)malloc(n * sizeof(double));
    random_vec(b, n);

    spsolve_LDLT(ilu0, x, b, n);
    mat_vec_csr_sym(mat, x, y);
    cblas_daxpy(n, -1., b, 1, y, 1);
    double err = cblas_dnrm2(n, y, 1) / cblas_dnrm2(n, b, 1);
    printf("ilu0 err = %.3le\n", err);

    spsolve_LDLT(ilu1, x, b, n);
    mat_vec_csr_sym(mat, x, y);
    cblas_daxpy(n, -1., b, 1, y, 1);
    err = cblas_dnrm2(n, y, 1) / cblas_dnrm2(n, b, 1);
    printf("ilu1 err = %.3le\n", err);

    free_csr(mat);
    free_csr(ilu0);
    free_csr(ilu1);
    free(b);
    free(x);
    free(y);
}

int main() {

    srand(10);
    test_csr_solve(12); // working
    test_ilu(10, 8, 6);

    return 0;
}