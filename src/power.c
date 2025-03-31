#include "power.h"
#include "matrix.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(a) ((a) * (a))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define PRECISION 10

/**
 * @brief Remplit un vecteur de taille n avec des valeurs aléatoires
 * @note Bonne pratique pour calculer un vecteur propre
 */
void randomize_eigv(double *v, int n) {
    for (int i = 0; i < n; i++) {
        v[i] = drand48();
        // v[i] = 1.;
    }
    double norm = cblas_dnrm2(n, v, 1);
    cblas_dscal(n, 1. / norm, v, 1);
}

/**
 * @brief Computes the lowest energy eigv/eigw of K x = lb M x
 * @param K Stiffness matrix in band storage, factorized LDL'
 * @param M Mass matrix in csr storage
 * @param eigw Pointer to the shift value, set to the lowest eigenvalue
 * @param eigv Pointer to the eigenvector (size n)
 * @param rel_tol Tolerance for the convergence
 * @param max_iter Maximum number of iterations
 * @param deflate Number of eigenvectors to deflate against
 * @param verbose 1 to print the convergence information, 0 otherwise
 * @return the number of iterations if converged, -1 otherwise
 */
int generalized_power(
    const SymBandMatrix *K,
    const CSRMatrix *M,
    double *eigw,
    double *eigv,
    const double rel_tol,
    const int max_iter,
    const int deflate,
    const int verbose
) {

    if (verbose) {
        printf("\n===========   Generalized Power iteration   ===========\n");
        printf("                shift = %11.5le\n", eigw[0]);
    }

    int iter, i;
    int k = K->k;
    int n = K->n;
    double diff_eigw, diff_eigv, lb, lb_prev, norm, eps;
    double *x = eigv;
    double *y = (double *)malloc(n * sizeof(double)); // previous
    double *z = (double *)malloc(n * sizeof(double)); // orthog. buffer
    double *LDL = K->data;

    randomize_eigv(x, n); // Initialize the eigenvector
    diff_eigw = 1.;       // Difference in eigenvalue
    diff_eigv = 1.;       // Difference in eigenvector
    lb = lb_prev = 0.;    // Previous eigenvalue guess
    eps = rel_tol;        // Convergence tolerance

    csr_sym_mat_vec(n, M, x, z);

    for (iter = 0; (eps < diff_eigw) && (iter < max_iter); iter++) {

        // Solve the linear system LDL' v = M v_prev
        cblas_dcopy(n, z, 1, x, 1);
        cblas_dcopy(n, z, 1, y, 1);
        solve_sym_band(LDL, n, k, x);
        
        // Orthogonalize (wrt M) the eigenvector if necessary
        for (i = 1; i <= deflate; i++) {
            csr_sym_mat_vec(n, M, x, z);
            norm = -cblas_ddot(n, eigv - i * n, 1, z, 1);
            cblas_daxpy(n, norm, eigv - i * n, 1, x, 1);
        }

        csr_sym_mat_vec(n, M, x, z);
        norm = sqrt(cblas_ddot(n, x, 1, z, 1)); // x' M x
        cblas_dscal(n, 1. / norm, x, 1);        // Normalize the eigenvector
        cblas_dscal(n, 1. / norm, z, 1);        // Normalize the eigenvector
        diff_eigv = cblas_ddot(n, x, 1, y, 1);  // Succ. eigvs M-orthogonality
        lb = eigw[0] + diff_eigv / norm;        // Rayleigh quotient
        diff_eigv = 1. - fabs(diff_eigv);       // Eigenvector

        eps = fmax(rel_tol * lb, 1e-14);
        diff_eigw = fabs(lb - lb_prev); // Eigenvalue convergence
        lb_prev = lb;                   // Update the previous eigenvalue

        if (verbose) {
            printf("iter: %3d | ", iter);
            printf("\033[1mλ = %15.10le\033[0m | ", lb);
            printf("Δλ = %9.3le  | ", diff_eigw);
            printf("Δv = %9.3le", diff_eigv);
            printf("\n");
        }
    }

    eigw[0] = lb;
    free(y);
    free(z);
    return (iter < max_iter) ? iter : -1;
}

double check_eig(const SymBandMatrix *K, CSRMatrix *M, double *x, double l) {
    int n = K->n;
    int k = K->k;
    double tmp;
    double *y = (double *)malloc(n * sizeof(double));

    if (M)
        csr_sym_mat_vec(n, M, x, y);
    else
        cblas_dcopy(n, x, 1, y, 1);
    tmp = cblas_ddot(n, x, 1, y, 1);

    cblas_dsbmv(BRM, BLW, n, k, 1., K->data, k + 1, x, 1, 0., y, 1);
    tmp = cblas_ddot(n, x, 1, y, 1) / tmp;

    free(y);
    return fabs(tmp - l);
}

/**
 * @brief Computes eigv/eigw near mu of a Generalized SDP Eigenproblem
 * @param K Stiffness matrix
 * @param M Mass matrix
 * @param eigw Eigenvalue pointer (size 1) that contains the shift mu as input
 * @param eigv Eigenvector array (ndofs)
 * @param rtol Relative tolerance for the convergence
 * @param max_it Maximum number of iterations
 */
void compute_eigv_shift(
    SymBandMatrix *K,
    SymBandMatrix *M,
    double *eigw,
    double *eigv,
    double rtol,
    int max_it
) {
    int n = K->n;
    int k = K->k;
    int ldK = k + 1;
    double mu = eigw[0];

    if (M) { // Compute A = K - mu M
        for (int i = 0; i < n * ldK; i++) {
            K->data[i] -= eigw[0] * M->data[i];
        }
    } else { // Compute A = K - mu I
        for (int i = 0; i < n; i++) {
            K->data[i * ldK + k] -= eigw[0];
        }
    }

    // Factorize the matrix in-place (in K)
    sym_band_LDL(K->data, n, k);

    // Transform M into CSR format
    CSRMatrix *M_csr = (M) ? band_to_csr(M) : NULL;

    // Compute the eigenvector
    int n_it = generalized_power(K, M_csr, eigw, eigv, rtol, max_it, 0, 0);
    if (M_csr)
        free_csr(M_csr);

    printf("μ = %*.*le | ", 10, 3, mu);
    if (n_it < 0) {
        printf("\033[1;31mλ = ");
        printf("%*.*le\033[0m | ", PRECISION + 7, PRECISION, eigw[0]);
        printf("did not converge in %d it !\n\n", max_it);
    } else {
        printf("\033[1mλ = %*.*le\033[0m", PRECISION + 7, PRECISION, eigw[0]);
        printf(" | converged in %4d it\n\n", n_it);
    }
}

/**
 * @brief Computes eigws and eigvs of a Generalized SDP Eigenproblem
 * @param K Stiffness matrix
 * @param M Mass matrix
 * @param eigws Eigenvalues array (nb)
 * @param eigvs Eigenvector array (nb x ndofs)
 * @param nb Number of eigenmodes to compute
 * @param rtol Relative tolerance for the convergence
 * @param max_it Maximum number of iterations
 */
void compute_eigvs_deflation(
    SymBandMatrix *K,
    SymBandMatrix *M,
    double *eigws,
    double *eigvs,
    int nb,
    double rtol,
    int max_it
) {

    int n = K->n;
    int k = K->k;
    int n_it;
    double *eigw, *eigv, err;
    SymBandMatrix *K_zero = allocate_sym_band_matrix(n, k);
    memcpy(K_zero->data, K->data, n * (k + 1) * sizeof(double));

    // Factorise la matrice K en LDL'
    double shift = 0.; // avoid zero eigws that make K singular
    for (int i = 0; i < n * (k + 1); i++) {
        K->data[i] -= shift * M->data[i];
    }
    sym_band_LDL(K->data, n, k);

    // Transforme M en format CSR
    CSRMatrix *M_csr = (M) ? band_to_csr(M) : NULL;

    for (int i = 0; i < nb; i++) {
        eigws[i] = shift;
        eigw = &eigws[i];
        eigv = &eigvs[i * n];

        n_it = generalized_power(K, M_csr, eigw, eigv, rtol, max_it, i, 0);
        if (n_it < 0) {
            printf("\033[1;31mλ[%2d] = ", i + 1);
            printf("%*.*le\033[0m | ", PRECISION + 7, PRECISION, eigws[i]);
            printf("did not converge in %d it !\n", max_it);
        } else {
            err = check_eig(K_zero, M_csr, eigv, eigws[i]);
            printf("\033[1mλ[%2d] = ", i + 1);
            printf("%*.*le\033[0m | ", PRECISION + 7, PRECISION, eigws[i]);
            printf("converged in %4d it | ", n_it);
            printf("error %9.3le\n", err);
        }
    }

    free_sym_band_matrix(K_zero);
    if (M_csr)
        free_csr(M_csr);
}

// -----------------------------------------------------------------------------
// TESTS
// -----------------------------------------------------------------------------

void lapack_soluce(
    const SymBandMatrix *K, const SymBandMatrix *M, double *eigws
) {
    int n = K->n;
    int k = K->k;

    double *K_copy = (double *)malloc(n * (k + 1) * sizeof(double));
    double *M_copy = (double *)malloc(n * (k + 1) * sizeof(double));
    memcpy(K_copy, K->data, n * (k + 1) * sizeof(double));
    memcpy(M_copy, M->data, n * (k + 1) * sizeof(double));
    memset(eigws, 0, n * sizeof(double));

    // REQUIRES LAPACK LIBRARY
    // clang-format off
    // LAPACKE_dsbgv(
    //     LAPACK_COL_MAJOR, 'N', 'U', n, k, k, 
    //     K_copy, k+1, M_copy, k+1, eigws, NULL, n
    // );
    // clang-format on

    free(K_copy);
    free(M_copy);
}

/**
 * @brief Prints the X' K X matrix
 * @param A Matrix A in band storage row-major (n x n) with bandwidth k
 * @param X Matrix X in col-major (n x nb)
 */
void print_XAX(
    const double *A, const double *X, int n, int k, int nb, char *msg
) {
    nb = MIN(nb, 10);
    double *tmp = (double *)malloc(n * sizeof(double));
    double dot;
    printf("%s = \n", msg);
    for (int j = 0; j < nb; j++) {
        const double *v = X + j * n;
        if (A)
            cblas_dsbmv(BRM, BLW, n, k, 1., A, k + 1, v, 1, 0., tmp, 1);
        else
            cblas_dcopy(n, v, 1, tmp, 1);
        printf("  ");
        for (int i = 0; i < j; i++) {
            dot = cblas_ddot(n, tmp, 1, X + i * n, 1);
            printf("%9.2le  ", dot);
        }
        dot = cblas_ddot(n, tmp, 1, X + j * n, 1);
        printf("\033[1m%9.2le\033[0m  ", dot);
        printf("\n");
    }
    printf("\n");
    free(tmp);
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

    printf("\n===========  Deflation algorithm verification  ===========\n");
    printf("%18s %10s %15s\n", "Me", "LAPACK", "Difference");
    for (int i = 0; i < n; i++) {
        if (i < nb) {
            eps = fmax(1e-12, r_eps * fabs(eigws[i]));
            diff = fabs(eigws[i] - eigws_ref[i]);
            sprintf(
                info, "%18.3lf %10.3lf %15.6le", eigws[i], eigws_ref[i], diff
            );
        } else {
            eps = 1.;
            diff = 0.;
            sprintf(info, "%18s %10.3lf %15s", "", eigws_ref[i], "");
        }
        if (eps < diff) {
            sprintf(msg, "%s%s%s", red_color, info, dft_color);
            count++;
        } else {
            sprintf(msg, "%s", info);
        }
        printf("%s\n", msg);
    }
    if (count) {
        printf("%33s%d %10s\n", "", count, "errors");
    }
    printf("\n");
}

// int main() {
int test_power() {
    int n = 400;
    int k = 40;
    int nb = 10;

    srand48(0);
    SymBandMatrix *K = allocate_sym_band_matrix(n, k);
    SymBandMatrix *M = allocate_sym_band_matrix(n, k);
    for (int i = 0; i < n * (k + 1); i++) {
        K->data[i] = 2. * drand48() - 1.;
        M->data[i] = (2. * drand48() - 1.) * 0.;
    }
    for (int i = 0; i < n; i++) {
        K->data[i * (k + 1) + k] += 2 * k;
        K->data[i * (k + 1) + k / 2] += k;
        M->data[i * (k + 1) + k] += 1.;
    }
    int i1 = 3;
    for (int j = MAX(i1 - k, 0); j <= i1; j++)
        K->a[i1][j] = 0.;
    for (int j = i1; j < MIN(i1 + k + 1, n); j++)
        K->a[j][i1] = 0.;
    // print_sym_band(K);
    // exit(0);
    // print_sym_band(M);

    // Required arrays
    double *eigws = (double *)malloc((nb + 1) * n * sizeof(double));
    double *eigvs = eigws + n;

    // Test arrays
    double *eigws_lapack = (double *)malloc(n * sizeof(double));
    double *tmp = (double *)malloc(nb * sizeof(double));
    double *K_zero = (double *)malloc(n * (k + 1) * sizeof(double));
    double *M_zero = (double *)malloc(n * (k + 1) * sizeof(double));
    memcpy(K_zero, K->data, n * (k + 1) * sizeof(double));
    memcpy(M_zero, M->data, n * (k + 1) * sizeof(double));

    // TEST DEFLATION
    lapack_soluce(K, M, eigws_lapack);
    compute_eigvs_deflation(K, NULL, eigws, eigvs, nb, 1e-12, 5000);
    compare_sols(eigws, eigws_lapack, 1e-8, nb, nb);

    // TEST SHIFT
    // eigws[0] = 20.;
    // compute_eigv_shift(K, NULL, &eigws[0], eigvs, 1e-12, 1000);

    print_XAX(K_zero, eigvs, n, k, 1, "X' K X");
    print_XAX(NULL, eigvs, n, k, nb, "X' X");
    print_XAX(M_zero, eigvs, n, k, nb, "X' M X");

    free_sym_band_matrix(K);
    free_sym_band_matrix(M);
    free(eigws);

    free(eigws_lapack);
    free(tmp);
    free(K_zero);
    free(M_zero);
    return 0;
}
