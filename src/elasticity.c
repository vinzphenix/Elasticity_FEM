/**
 * File:    elasticity.c
 * Author:  Vincent Degrooff
 * Created: 2025
 *
 * Description:
 *   Matrix assembly and boundary conditions
 * 
 * Project:
 *   FEM Simulation Toolkit for Linear Elasticity
 */

#include "elasticity.h"
#include "model.h"
#include <cblas.h>
#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compute_shape_functions(
    ElementType e_type,
    size_t *nq,
    double w[4],
    double xi[4][2],
    double phi[4][4],
    double dph[4][4][2],
    const int boundary
) {
    double a = sqrt(1. / 3.), b;
    double x, y;
    // clang-format off
    // Set quadrature points
    if (boundary) {
        nq[0] = 2;
        b = (e_type == TRI) ? 0.5 : 1.0;
        w[0] = b; w[1] = b;
        x = (e_type == TRI) ? 0.5 : 0.;
        y = (e_type == TRI) ? 0. : -1.;
        b = (e_type == TRI) ? 0.5 : 1.;
        xi[0][0] = x - b * a; xi[0][1] = y;
        xi[1][0] = x + b * a; xi[1][1] = y;
    } else if (e_type == TRI) {
        nq[0] = 3;
        a = 1. / 6.;
        b = 2. / 3.;
        w[0] = a; w[1] = a; w[2] = a; w[3] = 0.;
        xi[0][0] = a ; xi[0][1] = a;
        xi[1][0] = b ; xi[1][1] = a;
        xi[2][0] = a ; xi[2][1] = b;
        xi[3][0] = 0.; xi[3][1] = 0.;
    } else if (e_type == QUAD) {
        // nq[0] = 1;
        // w[0] = 4.;
        // xi[0][0] = 0.; xi[0][1] = 0.;
        nq[0] = 4;
        w[0] = 1.; w[1] = 1.; w[2] = 1.; w[3] = 1.;
        xi[0][0] = -a; xi[0][1] = -a;
        xi[1][0] = +a; xi[1][1] = -a;
        xi[2][0] = +a; xi[2][1] = +a;
        xi[3][0] = -a; xi[3][1] = +a;
    }
    // Set shape functions and derivatives
    if (e_type == TRI) {
        for (int i = 0; i < *nq; i++) {
            phi[i][0] = 1. - xi[i][0] - xi[i][1];
            phi[i][1] = xi[i][0];
            phi[i][2] = xi[i][1];
            phi[i][3] = 0.;
            dph[i][0][0] = -1.; dph[i][0][1] = -1.;
            dph[i][1][0] = +1.; dph[i][1][1] = +0.;
            dph[i][2][0] = +0.; dph[i][2][1] = +1.;
            dph[i][3][0] = +0.; dph[i][3][1] = +0.;
        }
    } else if (e_type == QUAD) {
        for (int i = 0; i < *nq; i++) {
            x = xi[i][0];
            y = xi[i][1];
            phi[i][0] = 0.25 * (1. - x) * (1. - y);
            phi[i][1] = 0.25 * (1. + x) * (1. - y);
            phi[i][2] = 0.25 * (1. + x) * (1. + y);
            phi[i][3] = 0.25 * (1. - x) * (1. + y);
            dph[i][0][0] = -0.25 * (1. - y); dph[i][0][1] = -0.25 * (1. - x);
            dph[i][1][0] = +0.25 * (1. - y); dph[i][1][1] = -0.25 * (1. + x);
            dph[i][2][0] = +0.25 * (1. + y); dph[i][2][1] = +0.25 * (1. + x);
            dph[i][3][0] = -0.25 * (1. + y); dph[i][3][1] = +0.25 * (1. - x);
        }
    } else {
        printf("Unknown element type: %d\n", e_type);
        exit(EXIT_FAILURE);
    }
    // clang-format on
}

void set_hooke_matrix(double E, double nu, int m_type, double h[4][4]) {
    // Non typical Hooke matrix ordering (allows flexibility for axisymmetric)
    // relates (s_xx, s_yy, s_xy) to (e_xx, e_yy, e_xy) for plane stress/strain
    // relates (s_rr, s_zz, s_rz, s_tt) to (e_rr, e_zz, e_rz, e_tt) for axisym.
    double a;
    memset(&h[0][0], 0, 16 * sizeof(double));
    if (m_type == PLANE_STRESS) {
        a = E / (1. - nu * nu);
        h[0][0] = h[1][1] = a;          // volumetric diagonal
        h[0][1] = h[1][0] = a * nu;     // volumetric off-diagonal
        h[2][2] = E / (2. * (1. + nu)); // shear modulus
    } else if (m_type == PLANE_STRAIN) {
        a = E / ((1. + nu) * (1. - 2. * nu));
        h[0][0] = h[1][1] = a * (1. - nu);
        h[0][1] = h[1][0] = a * nu;
        h[2][2] = E / (2. * (1. + nu));
    } else if (m_type == AXISYMMETRIC) {
        a = E / ((1. + nu) * (1. - 2. * nu));
        h[0][0] = h[1][1] = h[3][3] = a * (1. - nu);
        h[0][1] = h[1][0] = a * nu;
        h[2][2] = E / (2. * (1. + nu));
        h[0][3] = h[1][3] = h[3][0] = h[3][1] = nu * a;
    } else {
        printf("Unknown model type: %d\n", m_type);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Set the strain basis matrix B
 * @param n_loc Number of local nodes
 * @param phi Shape functions (nloc)
 * @param dph Shape functions derivatives (nloc, 2)
 * @param B Strain basis matrix (nloc, 4)
 * @note Unset entries should be set to zero by the caller
 */
void set_strain_basis(
    const int n_loc,
    const double phi[4],
    const double dphi[4][2],
    const double inv_r,
    double B[8][4]
) {
    for (int i = 0; i < n_loc; i++) {
        B[2 * i + 0][0] = B[2 * i + 1][2] = dphi[i][0];
        B[2 * i + 0][2] = B[2 * i + 1][1] = dphi[i][1];
        B[2 * i + 0][3] = phi[i] * inv_r;
    }
}

void set_local_stiffness_matrix(
    const int n_loc,
    const double phi[4],
    const double dphi[4][2],
    const double H[4][4],
    const double r,
    const double inv_r,
    double K[8][8]
) {
    const int n_dof = 2 * n_loc;
    double BT[8][4] = {0};
    set_strain_basis(n_loc, phi, dphi, inv_r, BT);

    // Compute [eps(phi_i) H eps(phi_j)]
    double BT_H[8][4];
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < 4; j++) {
            BT_H[i][j] = 0.;
            for (int k = 0; k < 4; k++) {
                BT_H[i][j] += r * BT[i][k] * H[k][j];
            }
        }
    }
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < n_dof; j++) {
            K[i][j] = 0.;
            for (int k = 0; k < 4; k++) {
                K[i][j] += BT_H[i][k] * BT[j][k];
            }
        }
    }
}

void set_local_mass_matrix(
    const int n_loc, const double phi[4], const double r, double M[8][8]
) {
    for (int i = 0; i < n_loc; i++) {
        for (int j = 0; j < n_loc; j++) {
            M[2 * i + 0][2 * j + 0] = r * phi[i] * phi[j]; // u - u
            M[2 * i + 1][2 * j + 1] = r * phi[i] * phi[j]; // v - v
            M[2 * i + 0][2 * j + 1] = 0.;                  // u - v
            M[2 * i + 1][2 * j + 0] = 0.;                  // v - u
        }
    }
}

void set_geometry(
    const int n_loc,
    const double x_node[4][2],
    const double phi[4],
    const double dph[4][2],
    double x_ptr[2],
    double det[1],
    double dphi[4][2]
) {
    double dxidx[2][2] = {{0., 0.}, {0., 0.}};
    x_ptr[0] = x_ptr[1] = 0.;
    for (int j = 0; j < n_loc; j++) {
        x_ptr[0] += x_node[j][0] * phi[j];
        x_ptr[1] += x_node[j][1] * phi[j];
        dxidx[0][0] += x_node[j][0] * dph[j][0];
        dxidx[0][1] += x_node[j][0] * dph[j][1];
        dxidx[1][0] += x_node[j][1] * dph[j][0];
        dxidx[1][1] += x_node[j][1] * dph[j][1];
    }
    det[0] = dxidx[0][0] * dxidx[1][1] - dxidx[0][1] * dxidx[1][0];
    for (int j = 0; j < n_loc; j++) {
        dphi[j][0] = +dxidx[1][1] * dph[j][0] - dxidx[1][0] * dph[j][1];
        dphi[j][1] = -dxidx[0][1] * dph[j][0] + dxidx[0][0] * dph[j][1];
    }
}

/**
 * @brief Compute the mass and stiffness matrices of a FE mesh
 */
void assemble_system(FE_Model *model) {
    int n_elem = model->n_elem;
    int nl = model->n_local;
    const size_t *elem_nodes = model->elem_nodes;
    int n_node = model->n_node;
    const double *coords = model->coords;
    const size_t *idx_map = model->idx_map;
    double nu = model->nu;

    int max_diff = model->node_band;
    model->M = allocate_band_sym(2 * n_node, 2 * max_diff + 1);
    model->K = allocate_band_sym(2 * n_node, 2 * max_diff + 1);
    SymBandMatrix *M = model->M;
    SymBandMatrix *K = model->K;

    size_t nq;
    size_t num[4] = {0};
    double w[4], xi[4][2], phi[4][4], dph[4][4][2];
    double x_node[4][2], dphi[4][2], x_loc[2], det, scaleM, scaleK;
    double hooke[4][4], K_loc[8][8], M_loc[8][8];
    double r, ir;

    set_hooke_matrix(1., nu, model->m_type, hooke);
    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 0);

    for (size_t i = 0; i < n_elem; i++) {
        const size_t *local_nodes = &elem_nodes[nl * i];
        SET_ELEM_INFO(nl, local_nodes, idx_map, coords, x_node, num);
        for (size_t q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            r = (model->m_type == AXISYMMETRIC) ? x_loc[0] : 1.;
            ir = (model->m_type == AXISYMMETRIC) ? det / x_loc[0] : 0.;
            set_local_stiffness_matrix(nl, phi[q], dphi, hooke, r, ir, K_loc);
            set_local_mass_matrix(nl, phi[q], r, M_loc);
            scaleK = w[q] / det;
            scaleM = w[q] * det;
            // Global matrices are numbered [u1, v1, u2, v2, u3, v3, u4, v4]
            for (size_t j = 0; j < 2 * nl; j++) {
                size_t row = 2 * num[j / 2] + (j % 2);
                for (size_t k = 0; k < 2 * nl; k++) {
                    size_t col = 2 * num[k / 2] + (k % 2);
                    if (col <= row) { // only fill lower part
                        K->a[row][col] += K_loc[j][k] * scaleK;
                        M->a[row][col] += M_loc[j][k] * scaleM;
                    }
                }
            }
        }
    }
    return;
}

/**
 * @brief Add the bulk source term to the rhs
 * @param model Finite element model
 * @param rhs Right-hand side of the system
 * @param bulk_source Function pointer to the source term evaluation
 */
void add_bulk_source(FE_Model *model, double *rhs) {
    size_t nq;
    size_t num[4] = {0};
    double w[4], xi[4][2], phi[4][4], dph[4][4][2];
    double x_node[4][2], dphi[4][2], x_loc[2], det, r;
    double f[2];

    int n_elem = model->n_elem;
    int nl = model->n_local;
    const size_t *elem_nodes = model->elem_nodes;
    const double *coords = model->coords;
    const size_t *idx_map = model->idx_map;

    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 0);

    for (size_t i = 0; i < n_elem; i++) {
        const size_t *local_nodes = &elem_nodes[nl * i];
        SET_ELEM_INFO(nl, local_nodes, idx_map, coords, x_node, num);
        for (size_t q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            r = (model->m_type == AXISYMMETRIC) ? x_loc[0] : 1.;
            model->set_bk_source(model->rho, x_loc, f);
            f[0] *= model->L_ref / model->E;
            f[1] *= model->L_ref / model->E;
            for (size_t j = 0; j < nl; j++) {
                rhs[2 * num[j] + 0] += w[q] * det * phi[q][j] * f[0] * r;
                rhs[2 * num[j] + 1] += w[q] * det * phi[q][j] * f[1] * r;
            }
        }
    }
}

double set_direction(
    const char kind, const double xy1[2], const double xy2[2], double *dir
) {
    double dx = xy2[0] - xy1[0];
    double dy = xy2[1] - xy1[1];
    double det = hypot(dx, dy);
    if (kind == 'x') {
        dir[0] = 1.0;
        dir[1] = 0.0;
    } else if (kind == 'y') {
        dir[0] = 0.0;
        dir[1] = 1.0;
    } else if (kind == 'n') {
        dir[0] = +dy / det;
        dir[1] = -dx / det;
    } else if (kind == 't') {
        dir[0] = dx / det;
        dir[1] = dy / det;
    } else { // cannot happen, but avoids compiler warning
        dir[0] = 0.0;
        dir[1] = 0.0;
    }
    return det;
}

/**
 * @brief Add the Robin boundary condition (Neumann is a special case)
 * @param model Finite element model
 * @param entity Entity number
 * @param size Size of edges (2 * n_edges)
 * @param edges Nodes on edges of the boundary (2 * n_edges)
 * @param kind x, y, n, or t for x, y, normal, or tangent
 * @param rhs System rhs to be modified
 * @note The Robin bc : n*sigma*kind = -k*(u - u_ref) + f = -alpha u + beta
 */
void apply_force(
    const FE_Model *model,
    const size_t entity,
    const size_t size,
    const size_t *edges,
    const char kind,
    double *rhs
) {
    size_t j1, j2, num1, num2, num_min, num_max, nq;
    double r, ref_l, det, val, dir[2], xy[2], f[2], tuu, tuv, tvv;
    double w[4], xi[4][2], phi[4][4], dph[4][4][2];
    double L = model->L_ref;
    double *x = model->coords;
    SymBandMatrix *K = model->K;

    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 1);
    ref_l = cblas_dsum(nq, w, 1);

    for (int e = 0; 2 * e < size; e += 1) {
        j1 = edges[2 * e + 0] - 1;
        j2 = edges[2 * e + 1] - 1;
        det = set_direction(kind, &x[2 * j1 + 0], &x[2 * j2 + 0], dir) / ref_l;
        num1 = 2 * model->idx_map[j1];
        num2 = 2 * model->idx_map[j2];
        num_min = MIN(num1, num2);
        num_max = MAX(num1, num2);
        for (int q = 0; q < nq; q++) {
            xy[0] = (x[2 * j1 + 0] * phi[q][0] + x[2 * j2 + 0] * phi[q][1]) * L;
            xy[1] = (x[2 * j1 + 1] * phi[q][0] + x[2 * j2 + 1] * phi[q][1]) * L;
            model->set_bd_force(entity, kind, xy, f);
            f[0] *= model->L_ref / model->E; // alpha (dimensionless)
            f[1] *= 1. / model->E;           // beta (dimensionless)
            r = (model->m_type == AXISYMMETRIC) ? xy[0] / L : 1.;
            // Boundary stiffness matrix
            tuu = f[0] * dir[0] * dir[0];
            tuv = f[0] * dir[0] * dir[1];
            tvv = f[0] * dir[1] * dir[1];
            val = w[q] * phi[q][0] * phi[q][0] * r * det; // j1 j1
            K->a[num1 + 0][num1 + 0] += val * tuu;
            K->a[num1 + 1][num1 + 0] += val * tuv;
            K->a[num1 + 1][num1 + 1] += val * tvv;
            val = w[q] * phi[q][1] * phi[q][1] * r * det; // j2 j2
            K->a[num2 + 0][num2 + 0] += val * tuu;
            K->a[num2 + 1][num2 + 0] += val * tuv;
            K->a[num2 + 1][num2 + 1] += val * tvv;
            val = w[q] * phi[q][1] * phi[q][0] * r * det; // j2 j1
            K->a[num_max + 0][num_min + 0] += val * tuu;
            K->a[num_max + 0][num_min + 1] += val * tuv;
            K->a[num_max + 1][num_min + 0] += val * tuv;
            K->a[num_max + 1][num_min + 1] += val * tvv;
            // Boundary force vector
            val = w[q] * r * det * f[1]; // rhs
            rhs[num1 + 0] += phi[q][0] * val * dir[0];
            rhs[num1 + 1] += phi[q][0] * val * dir[1];
            rhs[num2 + 0] += phi[q][1] * val * dir[0];
            rhs[num2 + 1] += phi[q][1] * val * dir[1];
        }
    }
}

void compute_node_normals(
    const size_t n_edge, const double *x, double *nx, double *ny, int periodic
) {
    double dx, dy;
    for (int e = 0; e < n_edge; e += 1) { // edge normals
        dx = x[2 * (e + 1) + 0] - x[2 * (e + 0) + 0];
        dy = x[2 * (e + 1) + 1] - x[2 * (e + 0) + 1];
        nx[e] = +dy;
        ny[e] = -dx;
    }
    nx[n_edge] = nx[n_edge - 1];
    ny[n_edge] = ny[n_edge - 1];
    for (int e = n_edge - 1; 0 < e; e--) {
        nx[e] = 0.5 * (nx[e - 1] + nx[e]);
        ny[e] = 0.5 * (ny[e - 1] + ny[e]);
    }
    if (periodic) {
        nx[0] = 0.5 * (nx[0] + nx[n_edge]);
        ny[0] = 0.5 * (ny[0] + ny[n_edge]);
    }
    for (int e = 0; e <= n_edge; e += 1) {
        double norm = hypot(nx[e], ny[e]);
        nx[e] /= norm;
        ny[e] /= norm;
    }
}

void project_system(
    size_t i, SymBandMatrix *K, double *rhs, double nx, double ny
) {
    size_t jxy, i_bound;
    size_t k_node = K->k / 2; // k = 2 * k_node + 1
    size_t ix = 2 * i + 0;
    size_t iy = 2 * i + 1;
    double aix, aiy, axy;
    // Rotate rows
    i_bound = (i < k_node) ? 0 : 2 * (i - k_node);
    for (jxy = i_bound; jxy < ix; jxy++) {
        aix = K->a[ix][jxy];
        aiy = K->a[iy][jxy];
        K->a[ix][jxy] = nx * aix + ny * aiy;
        K->a[iy][jxy] = ny * aix - nx * aiy;
    }
    i_bound = MIN(2 * (i + k_node + 1), K->n);
    for (jxy = iy + 1; jxy < i_bound; jxy++) {
        aix = K->a[jxy][ix];
        aiy = K->a[jxy][iy];
        K->a[jxy][ix] = nx * aix + ny * aiy;
        K->a[jxy][iy] = ny * aix - nx * aiy;
    }
    // Central 2x2 block
    aix = K->a[ix][ix];
    aiy = K->a[iy][ix];
    axy = K->a[iy][iy];
    K->a[ix][ix] = nx * nx * aix + 2 * nx * ny * aiy + ny * ny * axy;
    K->a[iy][iy] = ny * ny * aix - 2 * nx * ny * aiy + nx * nx * axy;
    K->a[iy][ix] = nx * ny * aix + (ny * ny - nx * nx) * aiy - nx * ny * axy;
    // Rhs
    aix = rhs[ix];
    aiy = rhs[iy];
    rhs[ix] = nx * aix + ny * aiy;
    rhs[iy] = ny * aix - nx * aiy;
}

/**
 * @brief Apply Dirichlet boundary conditions
 * @param model Finite element model
 * @param entity Entity number
 * @param size Size of edges (2 * n_edges)
 * @param edges Nodes on edges of the boundary (2 * n_edges)
 * @param kind x, y, n, or t for x, y, normal, or tangent
 * @param rhs System rhs to be modified
 */
void apply_dirichlet(
    const FE_Model *model,
    const size_t entity,
    const size_t size,
    size_t *edges,
    const char kind,
    double *rhs
) {
    double val;
    size_t node, num, numxy, i, i_bound;
    const size_t *i_map = model->idx_map;
    SymBandMatrix *K = model->K;
    size_t n_bd_node = size / 2 + 1;
    int periodic = edges[0] == edges[size - 1];

    double *nx = malloc(n_bd_node * sizeof(double));
    double *ny = malloc(n_bd_node * sizeof(double));
    double *xy = malloc(2 * n_bd_node * sizeof(double));

    for (int j = 0; j < n_bd_node; j++) {
        edges[j] = edges[MIN(2 * j, size - 1)];
        node = edges[j] - 1;
        xy[2 * j + 0] = model->coords[2 * node + 0];
        xy[2 * j + 1] = model->coords[2 * node + 1];
    }
    compute_node_normals(n_bd_node - 1, xy, nx, ny, periodic);
    cblas_dscal(2 * n_bd_node, model->L_ref, xy, 1);  // scale for set_bd_disp
    n_bd_node = periodic ? n_bd_node - 1 : n_bd_node;

    for (int j = 0; j < n_bd_node; j++) {
        num = i_map[edges[j] - 1];
        model->set_bd_disp(entity, kind, &xy[2 * j], &val);
        val /= model->L_ref;
        if (kind == 'n' || kind == 't')
            project_system(num, K, rhs, nx[j], ny[j]);
        numxy = (kind == 'x' || kind == 'n') ? 2 * num + 0 : 2 * num + 1;
        // Zero out row associated to dof "numxy"
        i_bound = (numxy < K->k) ? 0 : numxy - K->k;
        for (i = i_bound; i < numxy; i++) {
            rhs[i] -= K->a[numxy][i] * val;
            K->a[numxy][i] = 0.;
        }
        // Zero out col associated to dof "numxy"
        i_bound = MIN(numxy + K->k + 1, K->n);
        for (i = numxy + 1; i < i_bound; i++) {
            rhs[i] -= K->a[i][numxy] * val;
            K->a[i][numxy] = 0.;
        }
        K->a[numxy][numxy] = 1.;
        rhs[numxy] = val;
        if (kind == 'n' || kind == 't')
            project_system(num, K, rhs, nx[j], ny[j]);
    }
    free(nx);
    free(ny);
    free(xy);
}

/**
 * @brief Enforce boundary conditions
 * @param model Finite element model
 * @param rhs Right-hand side of the system
 */
void enforce_bd_conditions(FE_Model *model, double *rhs) {
    int ierr;
    int *dt_phys, *ent;
    size_t dt_phys_n, entity_n;
    size_t *bd_nodes;
    size_t bd_n; // n_edge; (to be used for improvement)
    char *name, bkind;

    gmshModelGetPhysicalGroups(&dt_phys, &dt_phys_n, 1, &ierr);
    // printf("Physical groups: %zu\n", dt_phys_n);
    for (size_t i = 0; i < dt_phys_n; i += 2) {
        gmshModelGetPhysicalName(dt_phys[i], dt_phys[i + 1], &name, &ierr);
        // printf("  Physical group: %s\n", name);
        bkind = name[strlen(name) - 1];
        if (bkind != 'x' && bkind != 'y' && bkind != 'n' && bkind != 't') {
            exit(EXIT_FAILURE);
        }
        gmshModelGetEntitiesForPhysicalName(name, &ent, &entity_n, &ierr);
        for (size_t e = 0; e < entity_n; e += 2) {
            gmshModelMeshGetElementEdgeNodes(
                1, &bd_nodes, &bd_n, ent[e + 1], 1, 0, 1, &ierr
            );
            // size_t n_edge = bd_n / 2;
            // printf("    Entity %d : %zu nodes\n", ent[e + 1], n_edge + 1);
            if (strncmp(name, "fix", 3) == 0) {
                apply_dirichlet(model, ent[e + 1], bd_n, bd_nodes, bkind, rhs);
            } else if (strncmp(name, "force", 5) == 0) {
                apply_force(model, ent[e + 1], bd_n, bd_nodes, bkind, rhs);
            } else {
                printf("Unknown boundary condition: %s\n", name);
                exit(EXIT_FAILURE);
            }
            free(bd_nodes);
        }
        free(name);
        free(ent);
    }
    free(dt_phys);
}
