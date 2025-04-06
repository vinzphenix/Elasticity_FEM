#ifndef MATRIX_H
#define MATRIX_H

#include <stdbool.h>
#include <stddef.h>

#define RWMJ CblasRowMajor
#define LOWR CblasLower
#define UNIT CblasUnit
#define NNUN CblasNonUnit
#define TRNS CblasTrans
#define NTRS CblasNoTrans

typedef size_t idx_t;
#define FMT_I "%zu"
// typedef int idx_t;
// #define FMT_I "%d"

typedef enum SpSlvFlag {
    CsrTrans,
    CsrNoTrans,
    CsrDiag,
    CsrInvDiag,
    CsrUnit,
    CsrNonUnit
} SpSlvFlag;

typedef enum LinearSolver {
    Band,
    CG_NoPrec,
    CG_Jacobi,
    CG_SSOR,
    CG_ILU0,
    CG_ILU1
} LinearSolver;

typedef struct SymBandMatrix {
    idx_t n, k;   // dimension de la matrice et largeur de bande
    double *data; // 1D array [m*(k+1)] avec les Aij de la partie inférieure
    double **a;   // 1D array [m] de pointeurs vers chaque ligne -> a[i][j]
} SymBandMatrix;

typedef struct CSRMatrix {
    idx_t n, nnz;   // dimension de la matrice et nombre d'éléments non nuls
    idx_t *row_ptr; // 1D array [n+1] : indices de début de chaque ligne
    idx_t *col_idx; // 1D array [nnz] : colonnes des nnz
    double *data;   // 1D array [nnz] : valeurs des nnz
} CSRMatrix;

typedef struct CSCMatrix {
    idx_t n, nnz;   // dimension de la matrice et nombre d'éléments non nuls
    idx_t *col_ptr; // 1D array [n+1] : indices de début de chaque colonne
    idx_t *row_idx; // 1D array [nnz] : lignes des nnz
    idx_t *idx_csr; // 1D array [nnz] : mapping vers la matrice CSR
} CSCMatrix;

typedef void (*PrecSolveFn)(const void *M, double *z, const double *r, idx_t n);

void print_zu_vector(idx_t n, idx_t *vec);
void print_vector(idx_t n, double *vec);
void print_vector_row(idx_t n, double *vec);
void print_band_sym(SymBandMatrix *M);
void print_csr(CSRMatrix *M);

void write_vector(const char *filename, double *v, idx_t n);
void write_band_sym(SymBandMatrix *M, double *rhs, const char *filename);
void write_csr(CSRMatrix *csr, double *rhs, const char *filename);

SymBandMatrix *allocate_band_sym(const idx_t n, const idx_t k);
void free_band_sym(SymBandMatrix *mat);
void fctrz_band_sym(SymBandMatrix *A);
void solve_band_sym(const SymBandMatrix *LDL, double *x);

CSRMatrix *allocate_csr_matrix(const idx_t n, const idx_t nnz);
CSCMatrix *csr_to_csc(const CSRMatrix *csr);
CSRMatrix *band_to_csr_sym(const SymBandMatrix *band);
void free_csr(CSRMatrix *csr);
void solve_csr(const CSRMatrix *A, double *x, SpSlvFlag mode, SpSlvFlag unit);
void mat_vec_csr(const CSRMatrix *A, const double *x, double *y);
void mat_vec_csr_sym(const CSRMatrix *csr, const double *x, double *y);

CSRMatrix *ilu0_symmetric(const CSRMatrix *A);
CSRMatrix *ilu1_symmetric(const CSRMatrix *csr);

void spsolve_LDLT(const void *A, double *z, const double *r, idx_t n);
void spsolve_SSOR(const void *A, double *z, const double *r, idx_t n);
void spsolve_Jacobi(const void *M, double *z, const double *r, idx_t n);
void spsolve_none(const void *M, double *z, const double *r, idx_t n);
const char *solver_name(enum LinearSolver solver);

int PCG(
    const CSRMatrix *A,
    const void *M,
    PrecSolveFn p,
    const double *b,
    double *x,
    double eps,
    int max_it
);
int solve_system(
    SymBandMatrix *A, LinearSolver prec, const double *b, double *x
);

#endif
