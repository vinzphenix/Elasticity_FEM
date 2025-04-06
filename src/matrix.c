#include "matrix.h"
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_PRINT 30
#define N_DECIMALS 3

void print_vector_row(idx_t n, double *vec) {
    if (n > 31)
        return;
    for (idx_t i = 0; i < n; i++) {
        printf("%10.2le ", vec[i]);
    }
    printf("\n");
}

void print_vector(idx_t n, double *vec) {
    for (idx_t i = 0; i < n; i++) {
        printf("%4zu : %22.15le\n", i, vec[i]);
    }
    printf("\n");
}

void print_zu_vector(idx_t n, idx_t *vec) {
    for (idx_t i = 0; i < n; i++) {
        printf("%4zu : %6zu\n", i, vec[i]);
    }
    printf("\n");
}

void print_band_sym(SymBandMatrix *M) {
    if (M->n > MAX_PRINT)
        return;
    printf("\nSymmetric band matrix\n");
    idx_t j_min;
    double eps = 1.e-12;
    idx_t k = M->k;
    for (idx_t i = 0; i < M->n; i++) {
        j_min = (i < k) ? 0 : i - k;
        for (idx_t j = 0; j < j_min; j++) {
            printf("%10s ", "");
        }
        for (idx_t j = j_min; j <= i; j++) {
            if (fabs(M->a[i][j]) < eps)
                printf("%10s ", ".");
            else
                printf("%10.2le ", M->a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_csr(CSRMatrix *M) {
    if (M->n > MAX_PRINT)
        return;
    printf("\nCSR matrix\n%7s", "");

    for (idx_t j = 0; j < M->n; j++)
        printf("%*zu ", N_DECIMALS + 7, j);
    printf("\n");
    for (idx_t i = 0; i < M->n; i++) {
        printf("%4zu : ", i);
        for (idx_t j = 0, k = M->row_ptr[i]; j <= i; j++) {
            if (j == M->col_idx[k]) {
                printf("%*.*le ", N_DECIMALS + 7, N_DECIMALS, M->data[k]);
                k++;
            } else {
                printf("%*s ", N_DECIMALS + 7, "");
            }
        }
        printf("\n");
    }
}

void write_vector(const char *filename, double *v, idx_t n) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    for (idx_t i = 0; i < n; i++) {
        fprintf(f, "%22.15le\n", v[i]);
    }
    fclose(f);
}

void write_band_sym(SymBandMatrix *M, double *rhs, const char *filename) {
    int save_rhs = (rhs != NULL);
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    fprintf(f, "# SYM_BAND " FMT_I " " FMT_I " %d\n", M->n, M->k, save_rhs);
    idx_t k = M->k;
    for (idx_t i = 0; i < M->n; i++) {
        for (idx_t j = 0; j < k + 1; j++) {
            fprintf(f, "%22.15le ", M->data[i * (k + 1) + j]);
        }
        if (save_rhs)
            fprintf(f, "%22.15le", rhs[i]);
        fprintf(f, "\n");
    }
    fclose(f);
}

void write_csr(CSRMatrix *csr, double *rhs, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    int save_rhs = (rhs != NULL);
    fprintf(f, "# CSR " FMT_I " " FMT_I " %d\n", csr->n, csr->nnz, save_rhs);
    for (idx_t i = 0; i < csr->n; i++) {
        fprintf(f, FMT_I " ", csr->row_ptr[i]);
        if (save_rhs)
            fprintf(f, "%22.15le", rhs[i]);
        fprintf(f, "\n");
    }
    for (idx_t i = 0; i < csr->nnz; i++) {
        fprintf(f, FMT_I " %22.15le\n", csr->col_idx[i], csr->data[i]);
    }
    fclose(f);
}

/**
 * @brief Alloue une matrice symétrique sous format bande compressé
 * @param n Taille de la matrice
 * @param k Largeur de bande
 * @return Pointeur vers la matrice
 * @note Les éléments Aij sont stockés dans un tableau 1D de taille n*(k+1)
 * @note Une liste de pointeurs 'a' permet d'accéder à A_ij avec a[i][j]
 */
SymBandMatrix *allocate_band_sym(const idx_t n, const idx_t k) {
    // To make things simple, we allocate k+1 entries for each row
    // (even though we need less than that for the k first and last rows).
    // We want to have a[i][i] == data[k + i*(k+1)]

    SymBandMatrix *mat = (SymBandMatrix *)malloc(sizeof(SymBandMatrix));
    mat->n = n;
    mat->k = k;
    mat->data = (double *)calloc(n * (k + 1), sizeof(double));
    mat->a = (double **)malloc(n * sizeof(double *));
    // This is the tricky part :-)
    // We want a[i][i] == data[(k+1)*i + k]
    // which means a[i] + i == data + (k+1)*i + k
    // and therefore a[i] == data + k + k*i
    for (idx_t i = 0; i < n; i++) {
        mat->a[i] = mat->data + k + k * i;
    }

    /*
    x : matrix off diagonal element
    = : matrix diagonal element
    + : memory allocated but not used
    Example for k = 3:
    + + + =
      + + x =
        + x x =
          x x x =
            x x x =
              x x x =
                x x x =
                  x x x =
                    x x x =
    */
    return mat;
}

void free_band_sym(SymBandMatrix *mat) {
    free(mat->a);
    free(mat->data);
    free(mat);
}

/**
 * @brief Alloue une matrice CSR
 * @param n Taille de la matrice
 * @param nnz Nombre d'éléments non nuls
 * @return Pointeur vers la matrice
 */
CSRMatrix *allocate_csr_matrix(const idx_t n, const idx_t nnz) {
    CSRMatrix *csr = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    csr->n = n;
    csr->nnz = nnz;
    csr->row_ptr = (idx_t *)malloc((n + 1) * sizeof(idx_t));
    csr->col_idx = (idx_t *)malloc(nnz * sizeof(idx_t));
    csr->data = (double *)malloc(nnz * sizeof(double));
    return csr;
}

/**
 * @brief Converts a symmetric band matrix to CSR format
 * @param band Pointer to the band matrix
 * @return Pointer to the CSR matrix
 */
CSRMatrix *band_to_csr_sym(const SymBandMatrix *band) {
    double val, thresh = 1e-15;
    idx_t j_min, nnz = 0, k = band->k;

    for (idx_t i = 0; i < band->n; i++) {
        j_min = (i < k) ? 0 : i - k;
        for (idx_t j = j_min; j <= i; j++)
            if (thresh < fabs(band->a[i][j]))
                nnz++;
    }

    CSRMatrix *csr = allocate_csr_matrix(band->n, nnz);
    idx_t row_idx = 0, col_idx = 0;
    csr->row_ptr[0] = 0;

    for (idx_t i = 0; i < band->n; i++) {
        csr->row_ptr[row_idx + 1] = csr->row_ptr[row_idx];
        j_min = (i < k) ? 0 : i - k;
        for (idx_t j = j_min; j <= i; j++) {
            val = band->a[i][j];
            if (thresh < fabs(val)) {
                csr->row_ptr[row_idx + 1]++;
                csr->col_idx[col_idx] = j;
                csr->data[col_idx] = val;
                col_idx++;
            }
        }
        row_idx++;
    }

    return csr;
}

void free_csr(CSRMatrix *csr) {
    free(csr->row_ptr);
    free(csr->col_idx);
    free(csr->data);
    free(csr);
}

/**
 * @brief Alloue une matrice CSR
 * @param n Taille de la matrice
 * @param nnz Nombre d'éléments non nuls
 * @return Pointeur vers la matrice
 */
CSCMatrix *csr_to_csc(const CSRMatrix *csr) {
    idx_t n = csr->n;
    idx_t nnz = csr->nnz;

    CSCMatrix *csc = (CSCMatrix *)malloc(sizeof(CSCMatrix));
    csc->n = n;
    csc->nnz = nnz;

    csc->col_ptr = (idx_t *)calloc((n + 1), sizeof(idx_t));
    for (idx_t p = 0; p < nnz; p++)
        csc->col_ptr[csr->col_idx[p] + 1]++;

    for (idx_t i = 0; i < n; i++)
        csc->col_ptr[i + 1] += csc->col_ptr[i];

    csc->row_idx = (idx_t *)malloc(nnz * sizeof(idx_t));
    csc->idx_csr = (idx_t *)malloc(nnz * sizeof(idx_t));
    for (idx_t i = 0; i < n; i++) {
        for (idx_t p = csr->row_ptr[i]; p < csr->row_ptr[i + 1]; p++) {
            idx_t j = csr->col_idx[p];
            idx_t q = csc->col_ptr[j];
            csc->row_idx[q] = i;
            csc->idx_csr[q] = p;
            csc->col_ptr[j]++;
        }
    }
    for (idx_t i = n - 1; i > 0; i--)
        csc->col_ptr[i] = csc->col_ptr[i - 1];
    csc->col_ptr[0] = 0;

    return csc;
}

void free_csc(CSCMatrix *csc) {
    free(csc->col_ptr);
    free(csc->row_idx);
    free(csc->idx_csr);
    free(csc);
}

/**
 * @brief LDL' in-place d'une matrice symétrique bande
 * @param matrix Matrice symétrique bande
 */
void fctrz_band_sym(SymBandMatrix *matrix) {
    double akk, coef;
    idx_t pivot; // Index of the pivot diagonal element
    idx_t idx_i; // Index of the element ... rows below the pivot
    idx_t idx_j; // Index of the element ... columns to the right of the pivot,
                 // but actually below due to symmetry
    idx_t i_max; // Number of remaining rows below the pivot
    idx_t lda;   // Stride of the band storage

    double *A = matrix->data;
    idx_t n = matrix->n;
    idx_t b = matrix->k;
    lda = b + 1;

    for (idx_t k = 0; k < n; k++) {
        pivot = k * lda + b;
        akk = A[pivot];
        i_max = (k + b) < n ? b : n - k - 1;
        for (idx_t i = 1; i <= i_max; i++) {
            idx_j = pivot + (lda - 1);
            idx_i = pivot + (lda - 1) * i;
            coef = A[idx_i] / akk;
            for (idx_t j = 1; j <= i; j++) {
                A[idx_i + j] -= coef * A[idx_j];
                idx_j += lda - 1;
            }
        }
        for (idx_t i = 1; i <= i_max; i++) {
            idx_i = pivot + (lda - 1) * i;
            A[idx_i] /= akk;
        }
    }
}

/**
 * @brief Résout le système symétrique L D L' x = b
 * @param LDL Contient L dans sa partie inférieure et D sur la diagonale
 * @param b Vecteur de droite
 * @param x Vecteur solution
 */
void solve_band_sym(const SymBandMatrix *M, double *x) {
    idx_t n = M->n;
    idx_t b = M->k;
    // idx_t bound, lda = b + 1;
    double *LDL = M->data;

    // L z = b (with/without blas)
    cblas_dtbsv(RWMJ, LOWR, NTRS, UNIT, n, b, LDL, b + 1, x, 1);
    // for (idx_t i = 0; i < n; i++) {
    //     bound = (i < b) ? b - i : 0;
    //     for (idx_t j = bound; j < b; j++)
    //         x[i] -= LDL[i * lda + j] * x[j + i - b];
    // }

    // D y = z
    for (idx_t i = 0; i < n; i++)
        x[i] /= LDL[i * (b + 1) + b];

    // L'x = y (with/without blas)
    cblas_dtbsv(RWMJ, LOWR, TRNS, UNIT, n, b, LDL, b + 1, x, 1);
    // for (idx_t ii = 0; ii < n; ii++) {
    //     idx_t i = n - 1 - ii;
    //     bound = (i < b) ? b - i : 0;
    //     for (idx_t j = bound; j < b; j++)
    //         x[j + i - b] -= LDL[i * lda + j] * x[i];
    // }
}

/**
 * @brief Produit matrice-vecteur pour une matrice CSR
 * @param A Matrice CSR
 * @param x Vecteur d'entrée
 * @param y Vecteur de sortie
 */
void mat_vec_csr(const CSRMatrix *A, const double *x, double *y) {
    for (idx_t i = 0; i < A->n; i++) {
        y[i] = 0.;
        for (idx_t p = A->row_ptr[i]; p < A->row_ptr[i + 1]; p++)
            y[i] += A->data[p] * x[A->col_idx[p]];
    }
}

/**
 * @brief Produit matrice-vecteur pour une matrice symétrique CSR
 * @param A Matrice CSR dont seule la partie inférieure est stockée
 * @param x Vecteur d'entrée
 * @param y Vecteur de sortie
 * @note Tous les éléments diagonaux sont supposés non-nuls
 */
void mat_vec_csr_sym(const CSRMatrix *A, const double *x, double *y) {
    idx_t row_start, row_end, i, j, k;
    for (i = 0; i < A->n; i++) {
        row_start = A->row_ptr[i];
        row_end = A->row_ptr[i + 1];
        y[i] = 0.;
        for (k = row_start; k < row_end - 1; k++) {
            j = A->col_idx[k];
            y[i] += A->data[k] * x[j];
            y[j] += A->data[k] * x[i];
        }
        y[i] += A->data[row_end - 1] * x[i];
    }
}

/**
 * @brief Résout le système triagulaire Lx = b ou Ux = b
 * @param A Matrice CSR
 * @param x Vecteur solution
 * @param mode Indique si la matrice est transposée / diagonale
 * @param unit Indique si la diagonale est unitaire
 * @note La matrice est supposée stockée triangulaire inférieure
 */
void solve_csr(const CSRMatrix *A, double *x, SpSlvFlag mode, SpSlvFlag unit) {
    const idx_t n = A->n;
    const idx_t *rows = A->row_ptr;
    const idx_t *cols = A->col_idx;
    const double *data = A->data;

    if (mode == CsrNoTrans) {
        if (unit == CsrUnit) {
            for (idx_t i = 1; i < n; i++)
                for (idx_t p = rows[i]; p < rows[i + 1] - 1; p++)
                    x[i] -= data[p] * x[cols[p]];
        } else if (unit == CsrNonUnit) {
            for (idx_t i = 0; i < n; i++) {
                for (idx_t p = rows[i]; p < rows[i + 1] - 1; p++)
                    x[i] -= data[p] * x[cols[p]];
                x[i] /= data[rows[i + 1] - 1];
            }
        }
    } else if (mode == CsrTrans) {
        if (unit == CsrUnit) {
            for (idx_t i = n - 1; i > 0; i--)
                for (idx_t p = rows[i + 1] - 2; p >= rows[i]; p--)
                    x[cols[p]] -= data[p] * x[i];
        } else if (unit == CsrNonUnit) {
            for (idx_t i = n - 1; i > 0; i--) {
                x[i] /= data[rows[i + 1] - 1];
                for (idx_t p = rows[i + 1] - 2; p >= rows[i]; p--)
                    x[cols[p]] -= data[p] * x[i];
            }
            x[0] /= data[0]; // outside loop because idx_t unsigned
        }
    } else if (mode == CsrDiag) {
        for (idx_t i = 0; i < A->n; i++)
            x[i] /= A->data[A->row_ptr[i + 1] - 1];
    } else if (mode == CsrInvDiag) {
        for (idx_t i = 0; i < A->n; i++)
            x[i] *= A->data[A->row_ptr[i + 1] - 1];
    }
}

static inline double get_aij(
    const idx_t *rows, const idx_t *cols, const double *L, idx_t i, idx_t j
) {
    // Looks bad, but FEM matrix have few non-zero elements on each row
    for (idx_t k = rows[i]; k < rows[i + 1]; k++) {
        if (cols[k] == j)
            return L[k];
        if (cols[k] > j)
            return 0.;
    }
    return 0.;
}

/**
 * @brief Factorisation incomplète une matrice creuse --- ILU(0)
 * @param A Matrice creuse
 */
void ilu0_csr_inplace(CSRMatrix *A) {
    const idx_t n = A->n;
    const idx_t *rows = A->row_ptr;
    const idx_t *cols = A->col_idx;
    double *L = A->data;
    double coef;

    for (idx_t i = 1; i < n; i++) {
        for (idx_t j = rows[i]; j < rows[i + 1] - 1; j++) {
            coef = L[j] / L[rows[cols[j] + 1] - 1];
            for (idx_t k = j + 1; k < rows[i + 1]; k++)
                L[k] -= coef * get_aij(rows, cols, L, cols[k], cols[j]);
        }
    }
    for (idx_t i = 1; i < n; i++)
        for (idx_t j = rows[i]; j < rows[i + 1] - 1; j++)
            L[j] /= L[rows[cols[j] + 1] - 1];
}

/**
 * @brief Factorisation incomplète une matrice creuse --- ILU(0)
 * @param A Matrice creuse
 * @return Matrice ILU(0)
 */
CSRMatrix *ilu0_symmetric(const CSRMatrix *A) {
    CSRMatrix *ilu = allocate_csr_matrix(A->n, A->nnz);
    memcpy(ilu->row_ptr, A->row_ptr, (A->n + 1) * sizeof(A->col_idx));
    memcpy(ilu->col_idx, A->col_idx, A->nnz * sizeof(A->col_idx));
    memcpy(ilu->data, A->data, A->nnz * sizeof(A->data));
    ilu0_csr_inplace(ilu);
    return ilu;
}

int compare_size_t(const void *a, const void *b) {
    return (*(idx_t *)a - *(idx_t *)b);
}

/**
 * @brief Factorisation incomplète une matrice creuse --- ILU(1)
 * @param A Matrice creuse
 * @return Matrice ILU(1)
 */
CSRMatrix *ilu1_symmetric(const CSRMatrix *csr) {

    // Pas satisfait de l'implementation peu élégante
    idx_t n = csr->n;
    idx_t *rows = csr->row_ptr;
    idx_t *cols = csr->col_idx;

    CSCMatrix *csc = csr_to_csc(csr);
    idx_t *rows_ = csc->row_idx;
    idx_t *cols_ = csc->col_ptr;
    // idx_t *idxs_ = csc->idx_csr; (unused)

    idx_t *ilu_rows = calloc(n + 1, sizeof(idx_t));
    idx_t *marker = calloc(n, sizeof(idx_t));

    // First pass: count the number of non-zero elements
    for (idx_t i = 0; i < n; i++) {
        for (idx_t jp = rows[i]; jp < rows[i + 1] - 1; jp++)
            marker[cols[jp]] = i;
        for (idx_t jp = rows[i]; jp < rows[i + 1] - 1; jp++) {
            idx_t j = cols[jp];
            for (idx_t kq = cols_[j] + 1; rows_[kq] < i; kq++) {
                idx_t k = rows_[kq];
                ilu_rows[i + 1] += (marker[k] != i);
                marker[k] = i;
            }
        }
        ilu_rows[i + 1] += ilu_rows[i] + (rows[i + 1] - rows[i]);
    }

    CSRMatrix *ilu = (CSRMatrix *)malloc(sizeof(CSRMatrix));
    ilu->n = n;
    ilu->nnz = ilu_rows[n];
    ilu->row_ptr = ilu_rows;
    ilu->col_idx = malloc((ilu->nnz + 1) * sizeof(idx_t)); // due to (**)
    ilu->data = calloc(ilu->nnz, sizeof(double));

    // Second pass : fill the ILU(1)
    for (idx_t i = 0; i < n; i++) {
        idx_t p = ilu->row_ptr[i];
        for (idx_t jp = rows[i]; jp < rows[i + 1]; jp++) {
            marker[cols[jp]] = i;
            ilu->col_idx[p++] = cols[jp];
        }
        for (idx_t jp = rows[i]; jp < rows[i + 1] - 1; jp++) {
            idx_t j = cols[jp];
            for (idx_t kq = cols_[j] + 1; rows_[kq] < i; kq++) {
                idx_t k = rows_[kq];
                ilu->col_idx[p] = k;   // write by default (**)
                p += (marker[k] != i); // move ptr only if new
                marker[k] = i;
            }
        }
        idx_t st = ilu->row_ptr[i];
        idx_t fn = ilu->row_ptr[i + 1];
        qsort(&ilu->col_idx[st], fn - st, sizeof(idx_t), compare_size_t);
        for (idx_t jp = rows[i], jq = st; jp < rows[i + 1]; jq++) {
            int match = (csr->col_idx[jp] == ilu->col_idx[jq]);
            ilu->data[jq] = csr->data[jp] * match;
            jp += match;
        }
    }

    free(marker);
    free_csc(csc);

    // Most of the function runtime
    ilu0_csr_inplace(ilu);
    return ilu;
}

/**
 * @brief Résout le système symétrique creux L D L' x = b
 * @param LDL Contient L dans sa partie inférieure et D sur la diagonale
 * @param b Vecteur de droite
 * @param x Vecteur solution
 */
void spsolve_LDLT(const void *A, double *z, const double *r, idx_t n) {
    CSRMatrix *LDLT = (CSRMatrix *)A;
    memcpy(z, r, n * sizeof(*z));
    solve_csr(LDLT, z, CsrNoTrans, CsrUnit);
    solve_csr(LDLT, z, CsrDiag, CsrNonUnit);
    solve_csr(LDLT, z, CsrTrans, CsrUnit);
}

void spsolve_SSOR(const void *A, double *z, const double *r, idx_t n) {
    CSRMatrix *LDLT = (CSRMatrix *)A;
    memcpy(z, r, n * sizeof(*z));
    solve_csr(LDLT, z, CsrNoTrans, CsrNonUnit);
    solve_csr(LDLT, z, CsrInvDiag, CsrNonUnit);
    solve_csr(LDLT, z, CsrTrans, CsrNonUnit);
}

void spsolve_Jacobi(const void *M, double *z, const double *r, idx_t n) {
    CSRMatrix *A = (CSRMatrix *)M;
    memcpy(z, r, n * sizeof(*z));
    solve_csr(A, z, CsrDiag, CsrNonUnit);
}

void spsolve_none(const void *M, double *z, const double *r, idx_t n) {
    memcpy(z, r, n * sizeof(*z));
}

int PCG(
    const CSRMatrix *A,
    const void *M,
    PrecSolveFn solve,
    const double *b,
    double *x,
    double eps,
    int max_it
) {

    int it = 0;
    idx_t n = A->n;
    double alpha, beta = 0;
    double *d = (double *)malloc(n * sizeof(double));
    double *Ad = (double *)malloc(n * sizeof(double));
    double *r = (double *)malloc(n * sizeof(double));
    double *z = (double *)malloc(n * sizeof(double));

    memset(x, 0, n * sizeof(*x)); // init x
    memcpy(r, b, n * sizeof(*r)); // r = b - Ax
    solve(M, z, r, n);            // solve Mz = r
    memcpy(d, z, n * sizeof(*d)); // d = z
    double norm_r_zero = cblas_dnrm2(n, r, 1);
    double rz = cblas_ddot(n, r, 1, z, 1);

    while (eps * norm_r_zero < cblas_dnrm2(n, r, 1)) {
        // printf("Iteration %4d : res = %.5e \n", it, cblas_dnrm2(n, r, 1));
        mat_vec_csr_sym(A, d, Ad);
        alpha = rz / cblas_ddot(n, d, 1, Ad, 1);
        cblas_daxpy(n, alpha, d, 1, x, 1);
        cblas_daxpy(n, -alpha, Ad, 1, r, 1);
        beta = 1. / rz;
        solve(M, z, r, n);
        rz = cblas_ddot(n, r, 1, z, 1);
        beta *= rz;
        cblas_dscal(n, beta, d, 1);
        cblas_daxpy(n, 1, z, 1, d, 1);
        it++;
    }

    // if (it < max_it) {
    //     rn = cblas_dnrm2(n, r, 1) / norm_2_ro;
    //     printf("PCG converged in %d it. : res = %.3le\n", it, rn);
    // } else {
    //     printf("\033[1;31mCG did not converge\033[0m in %d it.\n", it);
    // }

    free(d);
    free(Ad);
    free(r);
    free(z);
    return (it < max_it) ? it : -1;
}

const char *solver_name(enum LinearSolver solver) {
    switch (solver) {
    case Band:
        return "Band_LDLT";
    case CG_NoPrec:
        return "CG_NoPrec";
    case CG_Jacobi:
        return "CG_Jacobi";
    case CG_ILU0:
        return "CG_ILU0";
    case CG_ILU1:
        return "CG_ILU1";
    case CG_SSOR:
        return "CG_SSOR";
    default:
        return "Unknown";
    }
}

/**
 * @brief Résout le système symétrique Ax = b
 * @param A_band Matrice symétrique bande
 * @param b Vecteur de droite
 * @param solver LinearSolver type
 * @param x Vecteur solution
 * @return Number of iterations
 */
int solve_system(
    SymBandMatrix *A_band, LinearSolver solver, const double *b, double *x
) {

    if (solver == Band) {
        memcpy(x, b, A_band->n * sizeof(*x));
        fctrz_band_sym(A_band);
        solve_band_sym(A_band, x);
        return 1;
    } else {
        const CSRMatrix *A_csr = band_to_csr_sym(A_band);
        const CSRMatrix *M_csr = NULL;
        PrecSolveFn solve_fn = NULL;
        const double eps = 1e-10;
        const int max_it = 1e4;
        int it;

        if (solver == CG_NoPrec) {
            M_csr = NULL;
            solve_fn = spsolve_none;
        } else if (solver == CG_Jacobi) {
            M_csr = A_csr;
            solve_fn = spsolve_Jacobi;
        } else if (solver == CG_SSOR) {
            M_csr = A_csr;
            solve_fn = spsolve_SSOR;
        } else if (solver == CG_ILU0) {
            M_csr = ilu0_symmetric(A_csr);
            solve_fn = spsolve_LDLT;
        } else if (solver == CG_ILU1) {
            M_csr = ilu1_symmetric(A_csr);
            solve_fn = spsolve_LDLT;
        } else {
            fprintf(stderr, "Unknown solver type\n");
            return -1;
        }

        it = PCG(A_csr, (const void *)M_csr, solve_fn, b, x, eps, max_it);

        if (solver == CG_ILU0 || solver == CG_ILU1)
            free_csr((CSRMatrix *)M_csr);
        free_csr((CSRMatrix *)A_csr);

        return it;
    }
}
