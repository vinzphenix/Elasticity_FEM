#define _POSIX_C_SOURCE 199309L
#include "matrix.h"
#include <cblas.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TIMESP(t1, t2)                                                         \
    ((double)((t2).tv_sec - (t1).tv_sec) +                                     \
     1e-9 * ((double)((t2).tv_nsec - (t1).tv_nsec)))

int compare(const void *a, const void *b) {
    return (*(size_t *)a - *(size_t *)b);
}

size_t cube_mat_FDM(size_t nx, CSRMatrix **csr_ptr) {
    size_t n = nx * nx * nx;
    const size_t deg = 6;
    size_t nnz = n + deg * n / 2;
    CSRMatrix *mat = allocate_csr_matrix(n, nnz);
    size_t *rows = mat->row_ptr;
    size_t *cols = mat->col_idx;
    double *vals = mat->data;
    size_t ptr = 0;
    size_t col;
    size_t neighs[deg + 1] = {};
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < nx; j++) {
            for (size_t k = 0; k < nx; k++) {
                size_t idx = i * nx * nx + j * nx + k;
                neighs[0] = (i + nx - 1) % nx * nx * nx + j * nx + k;
                neighs[1] = i * nx * nx + (j + nx - 1) % nx * nx + k;
                neighs[2] = i * nx * nx + j * nx + (k + nx - 1) % nx;
                neighs[3] = i * nx * nx + j * nx + k;
                neighs[4] = i * nx * nx + j * nx + (k + nx + 1) % nx;
                neighs[5] = i * nx * nx + (j + nx + 1) % nx * nx + k;
                neighs[6] = (i + nx + 1) % nx * nx * nx + j * nx + k;
                qsort(neighs, deg + 1, sizeof(size_t), compare);
                rows[idx] = ptr;
                for (size_t d = 0; d < deg + 1; d++) {
                    if (neighs[d] == idx) {
                        cols[ptr] = idx;
                        vals[ptr++] = deg;
                    } else if (neighs[d] < idx) {
                        cols[ptr] = neighs[d];
                        vals[ptr++] = -1.;
                    }
                }
            }
        }
    }
    rows[n] = ptr;
    for (ptr = rows[n - 1]; ptr < rows[n] - 1; ptr++) {
        vals[ptr] = 0.;
    }
    *csr_ptr = mat;
    return n;
}

typedef struct {
    double v;
    size_t j;
} jv;

int compare_jv(const void *a, const void *b) {
    return ((*(jv *)a).j - (*(jv *)b).j);
}

double L_value(int dim, int dx, int dy, int dz) {
    if (dim == 3) {
        if (dx * dy * dz)
            return -1. / 6.;
        else if ((dx * dy) || (dy * dz) || (dz * dx))
            return -1. / 3.;
        else if (dx || dy || dz)
            return 0.;
        else
            return 16. / 3.;
    } else if (dim == 2) {
        if (dx * dy)
            return -1. / 3.;
        else if (dx || dy)
            return -1. / 3.;
        else
            return 8. / 3.;
    } else if (dim == 1) {
        if (dx)
            return -0.5;
        else
            return 1.0;
    }
}

size_t regular_grid_FEM(int dim, size_t nx, CSRMatrix **csr_ptr) {
    const size_t ny = (1 < dim) ? nx : 1;
    const size_t nz = (2 < dim) ? nx : 1;
    const size_t n = nx * ny * nz;
    const size_t deg = 3 * ((1 < dim) ? 3 : 1) * ((2 < dim) ? 3 : 1) - 1;
    const size_t nnz = n + deg * n / 2;
    CSRMatrix *mat = allocate_csr_matrix(n, nnz);
    size_t *rows = mat->row_ptr;
    size_t *cols = mat->col_idx;
    double *vals = mat->data;
    size_t ptr = 0;
    size_t col;
    int y_bound = (1 < dim) ? 1 : 0;
    int z_bound = (2 < dim) ? 1 : 0;
    jv *neighs = malloc((deg + 1) * sizeof(jv));
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                size_t idx = i * nz * ny + j * nz + k;
                size_t nn = 0;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -y_bound; dy <= y_bound; dy++) {
                        for (int dz = -z_bound; dz <= z_bound; dz++) {
                            neighs[nn].j = ((i + dx + nx) % nx) * nz * ny +
                                           ((j + dy + ny) % ny) * nz +
                                           ((k + dz + nz) % nz);
                            neighs[nn].v = L_value(dim, dx, dy, dz);
                            nn++;
                        }
                    }
                }
                qsort(neighs, deg + 1, sizeof(jv), compare_jv);
                rows[idx] = ptr;
                for (size_t d = 0; d < deg + 1; d++) {
                    if (neighs[d].j <= idx) {
                        cols[ptr] = neighs[d].j;
                        vals[ptr++] = neighs[d].v;
                    }
                }
            }
        }
    }
    rows[n] = ptr;
    for (ptr = rows[n - 1]; ptr < rows[n] - 1; ptr++) {
        vals[ptr] = 0.;
    }
    *csr_ptr = mat;
    free(neighs);
    return n;
}

void random_vec_grid(int dim, double *v, size_t nx) {
    const double pi = 3.14159265;
    const int ny = (1 < dim) ? nx : 1;
    const int nz = (2 < dim) ? nx : 1;
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                size_t idx = i * nz * ny + j * nz + k;
                // double arg = 3 * pi / nx;
                double xi = 2. * (double)i / (double)nx - 1.;
                double yi = (1 < dim) ? 2. * (double)j / (double)ny - 1. : 0.;
                double zi = (2 < dim) ? 2. * (double)k / (double)nz - 1. : 0.;
                double ri2 = xi * xi + yi * yi + zi * zi;
                v[idx] = 1e2 * exp(-ri2 / (2. * 0.2 * 0.2));
            }
        }
    }
}

int test_CG(CSRMatrix *mat, double *b, LinearSolver prec, char *filename) {

    int n_it;
    size_t n = mat->n;
    struct timespec t1, t2;
    void *M = NULL;
    PrecSolveFn solve;
    double dt1, dt2, err, extra;

    if (prec == CG_NoPrec) {
        dt1 = 0.;
        extra = 0.;
        solve = spsolve_none;
    } else if (prec == CG_Jacobi) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        M = (void *)(mat);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        solve = spsolve_Jacobi;
        dt1 = TIMESP(t1, t2);
        extra = 0.;
    } else if (prec == CG_SSOR) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        M = (void *)(mat);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        solve = spsolve_SSOR;
        dt1 = TIMESP(t1, t2);
        extra = 0.;
    } else if (prec == CG_ILU0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        M = (void *)ilu0_symmetric(mat);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        solve = spsolve_LDLT;
        dt1 = TIMESP(t1, t2);
        extra = 1e2 * ((CSRMatrix *)M)->nnz / (double)mat->nnz;
    } else if (prec == CG_ILU1) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        M = (void *)ilu1_symmetric(mat);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        solve = spsolve_LDLT;
        dt1 = TIMESP(t1, t2);
        extra = 1e2 * ((CSRMatrix *)M)->nnz / (double)mat->nnz;
    } else {
        fprintf(stderr, "Unknown preconditioner type\n");
        return -1;
    }

    double *x = malloc(mat->n * sizeof(double));
    double *y = malloc(mat->n * sizeof(double));

    clock_gettime(CLOCK_MONOTONIC, &t1);
    n_it = PCG(mat, M, solve, b, x, 1e-7, 1e5);
    clock_gettime(CLOCK_MONOTONIC, &t2);
    dt2 = TIMESP(t1, t2);

    // print_vector(mat->n, b);
    // print_vector(mat->n, x);

    mat_vec_csr_sym(mat, x, y);
    cblas_daxpy(mat->n, -1., b, 1, y, 1);
    err = cblas_dnrm2(mat->n, y, 1);
    err /= cblas_dnrm2(mat->n, b, 1);
    printf("%-9s |  ", solver_name(prec));
    printf("\033[1m%6.2lf s\033[0m  | ", dt1 + dt2);
    printf("res = %.3le | ", err);
    printf("%4d it.  | ", n_it);
    printf("%6.2lf s  | ", dt1);
    printf("memory + %6.2lf %%\n", extra);

    if (strcmp(filename, "") != 0) {
        FILE *f = fopen(filename, "a");
        if (f == NULL) {
            printf("Error opening file\n");
            return -1;
        }
        fprintf(f, "%7zu %8zu %9s ", mat->n, mat->nnz, solver_name(prec));
        fprintf(f, "%3d %.5le %.5le %.5le %.5le\n", n_it, err, dt1, dt2, extra);
        fclose(f);
    }

    free(x);
    free(y);
    if (prec == CG_ILU0 || prec == CG_ILU1) {
        free_csr(M);
    }
    return 0;
}

size_t n_to_nx(int dim, size_t n) {
    size_t nx;
    if (dim == 1) {
        nx = n;
    } else if (dim == 2) {
        nx = (size_t)sqrt(n);
    } else if (dim == 3) {
        nx = (size_t)cbrt(n);
    } else {
        fprintf(stderr, "Dimension must be between 1 and 3\n");
        nx = 0;
    }
    return nx;
}

void analyze_CG(int dim, int analyze) {

    size_t n, nx_min, nx_max;
    CSRMatrix *mat;
    double *b;
    double ratio;
    char filename[100] = "";

    if (dim == 1) {
        printf("Dimension : 1D\n");
        if (analyze)
            strcpy(filename, "newperf_CG_1D.txt");
    } else if (dim == 2) {
        printf("Dimension : 2D\n");
        if (analyze)
            strcpy(filename, "newperf_CG_2D.txt");
    } else if (dim == 3) {
        printf("Dimension : 3D\n");
        if (analyze)
            strcpy(filename, "newperf_CG_3D.txt");
    } else {
        fprintf(stderr, "Dimension must be between 1 and 3\n");
        return;
    }

    for (int run = 0; run < 6; run++) {

        // size_t nx = n_to_nx(dim, 1e3 * (1 << run)); // 1D (high kappa)
        size_t nx = n_to_nx(dim, 1e3 * (1 << run)); // 2D/3D
        n = regular_grid_FEM(dim, nx, &mat);
        // write_csr(mat, NULL, "new_matrix.txt");
        ratio = 2e2 * mat->nnz / (double)(n * (n + 1));
        print_csr(mat);

        b = malloc(mat->n * sizeof(double));
        random_vec_grid(dim, b, nx);

        printf("Matrix    | \033[1mn = %-6zu\033[0m |", mat->n);
        printf(" fill-in %.2lf %%  |", ratio);
        printf(" nnz = %.0le\n", (double)mat->nnz);

        test_CG(mat, b, CG_NoPrec, filename);
        test_CG(mat, b, CG_Jacobi, filename);
        test_CG(mat, b, CG_SSOR, filename);
        test_CG(mat, b, CG_ILU0, filename);
        test_CG(mat, b, CG_ILU1, filename);
        printf("\n--------------------------------------\n\n");

        free(b);
        free_csr(mat);
    }
}


// gcc ./perf_CG.c ../src/matrix.c -g -o exec -lopenblas -lm -fsanitize=address -I../include/
int main() {
    srand(0);

    // analyze_CG(1, 0);
    analyze_CG(2, 0);
    // analyze_CG(3, 0);

    return 0;
}
