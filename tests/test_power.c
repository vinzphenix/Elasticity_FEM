#include "lapacke.h"
#include "matrix.h"
#include "power.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void lapack_soluce(
    const SymBandMatrix *K, const SymBandMatrix *M, double *eigws
) {
    int n = K->n;
    int k = K->k;
    double *K_copy, *M_copy;

    memset(eigws, 0, n * sizeof(double));
    K_copy = (double *)malloc(n * (k + 1) * sizeof(double));
    memcpy(K_copy, K->data, n * (k + 1) * sizeof(double));

    // REQUIRES LAPACK LIBRARY
    // clang-format off
    if (M) {
        M_copy = (double *)malloc(n * (k + 1) * sizeof(double));
        memcpy(M_copy, M->data, n * (k + 1) * sizeof(double));
        LAPACKE_dsbgv(
            LAPACK_COL_MAJOR, 'N', 'U', n, k, k, 
            K_copy, k+1, M_copy, k+1, eigws, NULL, n
        );
        free(M_copy);
    } else {
        LAPACKE_dsbev(
            LAPACK_COL_MAJOR, 'N', 'U', n, k, 
            K_copy, k+1, eigws, NULL, n
        );
    }
    // clang-format on

    free(K_copy);
}

int cmp_abs(const void *a, const void *b) {
    double abs_a = fabs(*(double *)a);
    double abs_b = fabs(*(double *)b);
    return (abs_a > abs_b) - (abs_a < abs_b);
}

void compare_sols(
    double *eigws, double *eigws_ref, double r_eps, int nb, int n
) {
    double diff, eps;
    char red_color[] = "\033[1;31m";
    char dft_color[] = "\033[0m";
    char info[64], msg[128];
    int count = 0;
    qsort(eigws, nb, sizeof(double), cmp_abs);
    qsort(eigws_ref, n, sizeof(double), cmp_abs);

    printf("\n ----------  Deflation algorithm verification  ---------- \n");
    printf("|%17s %10s %15s %12s\n", "Me", "LAPACK", "Difference", "|");
    for (int i = 0; i < n; i++) {
        if (i < nb) {
            eps = fmax(1e-12, r_eps * fabs(eigws[i]));
            diff = fabs(eigws[i] - eigws_ref[i]);
            sprintf(
                info, "%17.3lf %10.3lf %15.6le", eigws[i], eigws_ref[i], diff
            );
        } else {
            eps = 1.;
            diff = 0.;
            sprintf(info, "%17s %10.3lf %15s", "", eigws_ref[i], "");
        }
        if (eps < diff) {
            sprintf(msg, "|%s%s%s %12s", red_color, info, dft_color, "|");
            count++;
        } else {
            sprintf(msg, "|%s %12s", info, "|");
        }
        printf("%s\n", msg);
    }
    if (count) {
        printf("|%32s%d %9s %12s\n", "", count, "errors", "|");
    }
    printf(" -------------------------------------------------------- \n");
}

void create_matrices(int n, int k, SymBandMatrix **Kp, SymBandMatrix **Mp) {
    srand48(12);
    SymBandMatrix *K = allocate_band_sym(n, k);
    SymBandMatrix *M = allocate_band_sym(n, k);
    for (int i = 0; i < n * (k + 1); i++) {
        K->data[i] = 2. * drand48() - 1.;
        M->data[i] = (2. * drand48() - 1.);
    }
    for (int i = 0; i < n; i++) {
        K->data[i * (k + 1) + k] += 2 * k;
        K->data[i * (k + 1) + k / 2] += k;
        M->data[i * (k + 1) + k] += k; // diagonal dominance
    }

    // Make K singular
    // int i1 = 3;
    // for (int j = MAX(i1 - k, 0); j <= i1; j++)
    //     K->a[i1][j] = 0.;
    // for (int j = i1; j < MIN(i1 + k + 1, n); j++)
    //     K->a[j][i1] = 0.;

    // print_band_sym(K);
    // print_band_sym(M);
    *Kp = K;
    *Mp = M;
}

int test_power(int n, int k, int nb, int M_flag, int info) {

    // Allocate matrices
    SymBandMatrix *K = NULL;
    SymBandMatrix *M = NULL;
    create_matrices(n, k, &K, &M);
    char *mode = (0 == nb) ? "Shift" : "Deflation";
    nb = MAX(nb, 1);

    // Required arrays
    double *eigws = (double *)malloc(nb * (n + 1) * sizeof(double));
    double *eigvs = eigws + nb;

    // Test arrays
    double *eigws_lapack = (double *)malloc(n * sizeof(double));
    double *tmp = (double *)malloc(nb * sizeof(double));
    double *K_init = (double *)malloc(n * (k + 1) * sizeof(double));
    double *M_init = (double *)malloc(n * (k + 1) * sizeof(double));
    memcpy(K_init, K->data, n * (k + 1) * sizeof(double));
    memcpy(M_init, M->data, n * (k + 1) * sizeof(double));

    SymBandMatrix *M_test = (M_flag) ? M : NULL;

    if (strcmp(mode, "Shift") == 0) {
        eigws[0] = 5.;
        compute_eigv_shift(K, M_test, &eigws[0], eigvs, 1e-12, 1000);
    } else if (strcmp(mode, "Deflation") == 0) {
        lapack_soluce(K, M_test, eigws_lapack);
        compute_eigvs_deflation(K, M_test, eigws, eigvs, nb, 1e-12, 5000);
        compare_sols(eigws, eigws_lapack, 1e-8, nb, nb);
    }

    if (info) {
        print_XAX(K_init, eigvs, n, k, nb, "X' K X");
        print_XAX(M_init, eigvs, n, k, nb, "X' X");
        print_XAX(M_init, eigvs, n, k, nb, "X' M X");
    }

    free_band_sym(K);
    free_band_sym(M);
    free(eigws);

    free(eigws_lapack);
    free(tmp);
    free(K_init);
    free(M_init);
    return 0;
}

int main() {
    int info = 0;

    printf("============================================================\n");
    printf("===================  Power method tests  ===================\n");
    printf("============================================================\n\n");

    printf("=================  EIGENPROBLEM K x = l x  =================\n\n");
    printf("-----------------   Inverse power shift   ------------------\n");
    test_power(20, 6, 0, 0, info);
    printf("-----------------        Deflation        ------------------\n");
    test_power(20, 6, 10, 0, info);
    printf("\n");

    printf("================  Eigenproblem K x = l M x  ================\n\n");
    printf("-----------------   Inverse power shift   ------------------\n");
    test_power(20, 6, 0, 1, info);
    printf("-----------------        Deflation        ------------------\n");
    test_power(20, 6, 10, 1, info);
    printf("\n");

    return 0;
}