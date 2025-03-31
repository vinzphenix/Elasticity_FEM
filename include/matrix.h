#ifndef MATRIX_H
#define MATRIX_H

#include <stdbool.h>
#include <stddef.h>

#define BRM CblasRowMajor
#define BLW CblasLower
#define BUN CblasUnit
#define BNU CblasNonUnit
#define BTR CblasTrans
#define BNT CblasNoTrans

typedef struct SymBandMatrix {
    size_t n, k;  // dimension de la matrice et largeur de bande
    double *data; // 1D array [m*(k+1)] avec les Aij de la partie inférieure
    double **a;   // 1D array [m] de pointeurs vers chaque ligne -> a[i][j]
} SymBandMatrix;

typedef struct CSRMatrix {
    size_t n, nnz;
    size_t *row_ptr;
    size_t *col_idx;
    double *data;
} CSRMatrix;

SymBandMatrix *allocate_sym_band_matrix(const size_t n, const size_t k);
void print_sym_band(SymBandMatrix *M);
void free_sym_band_matrix(SymBandMatrix *mat);
void write_sym_band(SymBandMatrix *M, double *rhs, const char *filename);
void write_csr(CSRMatrix *csr, double *rhs, const char *filename);
void write_vector(const char *filename, double *v, size_t n);
void print_vector_row(size_t n, double *vec);
void print_vector(size_t n, double *vec);

CSRMatrix *band_to_csr(const SymBandMatrix *band);
void csr_sym_mat_vec(
    size_t n, const CSRMatrix *csr, const double *x, double *y
);
void free_csr(CSRMatrix *csr);

void sym_band_LDL(double *A, size_t n, size_t b);
void solve_sym_band(double *L, size_t n, size_t b, double *x);

#endif