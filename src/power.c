/**
 * File:    power.c
 * Author:  Vincent Degrooff
 * Created: 2025
 *
 * Description:
 *   Power iteration for the modal analysis
 * 
 * Project:
 *   FEM Simulation Toolkit for Linear Elasticity
 */

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


void mat_vec_call(size_t n, const CSRMatrix *csr, const double *x, double *y) {
    if (!csr) {
        cblas_dcopy(n, x, 1, y, 1);
        return;
    } else {
        mat_vec_csr_sym(csr, x, y);
    }
}

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
    int n = K->n;
    double diff_eigw, diff_eigv, lb, lb_prev, norm, eps;
    double *x = eigv;
    double *y = (double *)malloc(n * sizeof(double)); // previous
    double *z = (double *)malloc(n * sizeof(double)); // orthog. buffer

    randomize_eigv(x, n); // Initialize the eigenvector
    diff_eigw = 1.;       // Difference in eigenvalue
    diff_eigv = 1.;       // Difference in eigenvector
    lb = lb_prev = 0.;    // Previous eigenvalue guess
    eps = rel_tol;        // Convergence tolerance

    mat_vec_call(n, M, x, z);

    for (iter = 0; (eps < diff_eigw) && (iter < max_iter); iter++) {

        // Solve the linear system LDL' v = M v_prev
        cblas_dcopy(n, z, 1, x, 1);
        cblas_dcopy(n, z, 1, y, 1);
        solve_band_sym(K, x);
        
        // Orthogonalize (wrt M) the eigenvector if necessary
        for (i = 1; i <= deflate; i++) {
            mat_vec_call(n, M, x, z);
            norm = -cblas_ddot(n, eigv - i * n, 1, z, 1);
            cblas_daxpy(n, norm, eigv - i * n, 1, x, 1);
        }

        mat_vec_call(n, M, x, z);
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

    mat_vec_call(n, M, x, y);
    tmp = cblas_ddot(n, x, 1, y, 1);

    cblas_dsbmv(RWMJ, LOWR, n, k, 1., K->data, k + 1, x, 1, 0., y, 1);
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
    fctrz_band_sym(K);

    // Transform M into CSR format
    CSRMatrix *M_csr = (M) ? band_to_csr_sym(M) : NULL;

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
    SymBandMatrix *K_zero = allocate_band_sym(n, k);
    memcpy(K_zero->data, K->data, n * (k + 1) * sizeof(double));

    // Factorise la matrice K en LDL'
    double shift = 0.; // avoid zero eigws that make K singular
    if (M) {
        cblas_daxpy(n * (k + 1), -shift, M->data, 1, K->data, 1);
    } else {
        for (int i = 0; i < n; i++) {
            K->data[i * (k + 1) + k] -= shift;
        }
    }
    
    fctrz_band_sym(K);

    // Transforme M en format CSR
    CSRMatrix *M_csr = (M) ? band_to_csr_sym(M) : NULL;

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

    free_band_sym(K_zero);
    if (M_csr)
        free_csr(M_csr);
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
            cblas_dsbmv(RWMJ, LOWR, n, k, 1., A, k + 1, v, 1, 0., tmp, 1);
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
