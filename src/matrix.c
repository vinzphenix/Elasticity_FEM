#include "matrix.h"
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

SymBandMatrix *allocate_sym_band_matrix(const size_t n, const size_t k) {
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
    for (size_t i = 0; i < n; i++) {
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

void print_sym_band(SymBandMatrix *M) {
    if (M->n > 31)
        return;
    printf("\nSymmetric band matrix\n");
    size_t j_min;
    double eps = 1.e-12;
    size_t k = M->k;
    for (size_t i = 0; i < M->n; i++) {
        j_min = (i < k) ? 0 : i - k;
        for (size_t j = 0; j < j_min; j++) {
            printf("%10s ", "");
        }
        for (size_t j = j_min; j <= i; j++) {
            if (fabs(M->a[i][j]) < eps)
                printf("%10s ", ".");
            else
                printf("%10.2le ", M->a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void write_sym_band(SymBandMatrix *M, double *rhs, const char *filename) {
    int save_rhs = (rhs != NULL);
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    fprintf(f, "# SYM_BAND %zu %zu %d\n", M->n, M->k, save_rhs);
    size_t k = M->k;
    for (size_t i = 0; i < M->n; i++) {
        for (size_t j = 0; j < k + 1; j++) {
            fprintf(f, "%22.15le ", M->data[i * (k + 1) + j]);
        }
        if (save_rhs)
            fprintf(f, "%22.15le", rhs[i]);
        fprintf(f, "\n");
    }
    fclose(f);
}

void free_sym_band_matrix(SymBandMatrix *mat) {
    free(mat->a);
    free(mat->data);
    free(mat);
}

void print_vector_row(size_t n, double *vec) {
    if (n > 31)
        return;
    for (size_t i = 0; i < n; i++) {
        printf("%10.2le ", vec[i]);
    }
    printf("\n");
}

void print_vector(size_t n, double *vec) {
    for (size_t i = 0; i < n; i++) {
        printf("%4zu : %22.15le\n", i, vec[i]);
    }
    printf("\n");
}

CSRMatrix *band_to_csr(const SymBandMatrix *band) {
    CSRMatrix *csr = malloc(sizeof(CSRMatrix));
    double thresh = 1e-15;
    double val;
    size_t j_min;
    size_t k = band->k;
    csr->n = band->n;
    csr->nnz = 0;

    for (size_t i = 0; i < band->n; i++) {
        j_min = (i < k) ? 0 : i - k;
        for (size_t j = j_min; j <= i; j++)
            if (thresh < fabs(band->a[i][j]))
                csr->nnz++;
    }

    size_t row_idx = 0;
    size_t col_idx = 0;
    csr->row_ptr = malloc((band->n + 1) * sizeof(size_t));
    csr->col_idx = malloc(csr->nnz * sizeof(size_t));
    csr->data = malloc(csr->nnz * sizeof(double));
    csr->row_ptr[0] = 0;

    for (size_t i = 0; i < band->n; i++) {
        csr->row_ptr[row_idx + 1] = csr->row_ptr[row_idx];
        j_min = (i < k) ? 0 : i - k;
        for (size_t j = j_min; j <= i; j++) {
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

void write_csr(CSRMatrix *csr, double *rhs, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    int save_rhs = (rhs != NULL);
    fprintf(f, "# CSR %zu %zu %d\n", csr->n, csr->nnz, save_rhs);
    for (size_t i = 0; i < csr->n; i++) {
        fprintf(f, "%zu ", csr->row_ptr[i]);
        if (save_rhs)
            fprintf(f, "%22.15le", rhs[i]);
        fprintf(f, "\n");
    }
    for (size_t i = 0; i < csr->nnz; i++) {
        fprintf(f, "%zu %22.15le\n", csr->col_idx[i], csr->data[i]);
    }
    fclose(f);
}

void write_vector(const char *filename, double *v, size_t n) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    for (size_t i = 0; i < n; i++) {
        fprintf(f, "%22.15le\n", v[i]);
    }
    fclose(f);
}

/**
 * @brief Produit matrice-vecteur pour une matrice symétrique stockée en CSR
 * @param csr Matrice stockée en CSR
 * @param x Vecteur d'entrée
 * @param y Vecteur de sortie
 * @note La matrice est supposée symétrique et avec diagonale non-nulle
 */
void csr_sym_mat_vec(
    size_t n, const CSRMatrix *csr, const double *x, double *y
) {
    if (!csr) {
        cblas_dcopy(n, x, 1, y, 1);
        return;
    }
    size_t row_start, row_end, i, j, k;
    for (i = 0; i < csr->n; i++) {
        row_start = csr->row_ptr[i];
        row_end = csr->row_ptr[i + 1];
        y[i] = 0.;
        for (k = row_start; k < row_end - 1; k++) {
            j = csr->col_idx[k];
            y[i] += csr->data[k] * x[j];
            y[j] += csr->data[k] * x[i];
        }
        // Assuming that diagonal elements are stored last and are non-zero
        y[i] += csr->data[row_end - 1] * x[i];
    }
}

/**
 * @brief LDL' in-place d'une matrice symétrique bande
 */
void sym_band_LDL(double *A, size_t n, size_t b) {
    double akk, coef;
    size_t pivot; // Index of the pivot diagonal element
    size_t idx_i; // Index of the element ... rows below the pivot
    size_t idx_j; // Index of the element ... columns to the right of the pivot,
                  // but actually below due to symmetry
    size_t i_max; // Number of remaining rows below the pivot
    size_t lda;   // Stride of the band storage

    lda = b + 1;
    for (size_t k = 0; k < n; k++) {
        pivot = k * lda + b;
        akk = A[pivot];
        i_max = (k + b) < n ? b : n - k - 1;
        for (size_t i = 1; i <= i_max; i++) {
            idx_j = pivot + (lda - 1);
            idx_i = pivot + (lda - 1) * i;
            coef = A[idx_i] / akk;
            for (size_t j = 1; j <= i; j++) {
                A[idx_i + j] -= coef * A[idx_j];
                idx_j += lda - 1;
            }
        }
        for (size_t i = 1; i <= i_max; i++) {
            idx_i = pivot + (lda - 1) * i;
            A[idx_i] /= akk;
        }
    }
}

void solve_sym_band(double *L, size_t n, size_t b, double *x) {
    // size_t bound, lda = b + 1;

    cblas_dtbsv(BRM, BLW, BNT, BUN, n, b, L, b + 1, x, 1);
    // for (size_t i = 0; i < n; i++) { 
    //     bound = (i < b) ? b - i : 0;
    //     for (size_t j = bound; j < b; j++) {
    //         x[i] -= L[i * lda + j] * x[j + i - b];
    //     }
    // }
    
    for (size_t i = 0; i < n; i++) {
        x[i] /= L[i * (b + 1) + b];
    }
    
    cblas_dtbsv(BRM, BLW, BTR, BUN, n, b, L, b + 1, x, 1);
    // for (size_t ii = 0; ii < n; ii++) {
    //     size_t i = n - 1 - ii;
    //     bound = (i < b) ? b - i : 0;
    //     for (size_t j = bound; j < b; j++) {
    //         x[j + i - b] -= L[i * lda + j] * x[i];
    //     }
    // }
}
