#include "visualize.h"
#include "elasticity.h"
#include "model.h"
#include "renumber.h"
#include <FL/math.h>
#include <cblas.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ZERO_RADIUS 1e-14
#define UNZIP6(a) &(a[0]), &(a[1]), &(a[2]), &(a[3]), &(a[4]), &(a[5])
#define PLOT gmshViewAddHomogeneousModelData
#define UPDATE_MIN(a, b) ((a) = (a) < (b) ? (a) : (b))
#define UPDATE_MAX(a, b) ((a) = (a) > (b) ? (a) : (b))

#define CART_TO_POLAR 0
#define DISPLAY_AVG 0

void compute_local_stress(
    const int n_loc,
    const double phi[4],
    const double dphi[4][2],
    const double u[8],
    const double h[4][4],
    const double ir,
    double stress[4]
) {
    double BT[8][4] = {0};
    double strain[4] = {0};
    set_strain_basis(n_loc, phi, dphi, ir, BT);
    for (int j = 0; j < 2 * n_loc; j++) {
        strain[0] += BT[j][0] * u[j];
        strain[1] += BT[j][1] * u[j];
        strain[2] += BT[j][2] * u[j];
        strain[3] += BT[j][3] * u[j];
    }
    for (int i = 0; i < 4; i++) {
        stress[i] = 0.;
        for (int j = 0; j < 4; j++) {
            stress[i] += h[i][j] * strain[j];
        }
    }
}

void cartesian_to_polar(size_t n_node, double *sigma, double *x) {
    double s_xx, s_yy, s_xy, r, c, s, c2, s2;
    for (size_t i = 0; i < n_node; i++) {
        s_xx = sigma[9 * i + 0];
        s_yy = sigma[9 * i + 4];
        s_xy = sigma[9 * i + 1];
        r = hypot(x[2 * i + 0], x[2 * i + 1]);
        c = x[2 * i + 0] / r;
        s = x[2 * i + 1] / r;
        c2 = c * c;
        s2 = s * s;
        sigma[9 * i + 0] = s_xx * c2 + s_yy * s2 + 2. * s_xy * c * s;
        sigma[9 * i + 4] = s_xx * s2 + s_yy * c2 - 2. * s_xy * c * s;
        sigma[9 * i + 1] = (s_yy - s_xx) * c * s + s_xy * (c2 - s2);
        sigma[9 * i + 3] = sigma[9 * i + 1];
    }
}

/**
 * @brief Compute the mass matrix for the stress least squares problem
 * @param model Finite element model
 * @param M Mass matrix with band storage (n_node x (k + 1))
 */
void assemble_mass_lsq(FE_Model *model, SymBandMatrix *M) {
    int n_elem = model->n_elem;
    int nl = model->n_local;
    const size_t *elem_nodes = model->elem_nodes;
    const double *coords = model->coords;
    const size_t *idx_map = model->idx_map;

    size_t nq;
    size_t num[4] = {0};
    double w[4], xi[4][2], phi[4][4], dph[4][4][2], r;
    double x_node[4][2], x_loc[2], dphi[4][2], det, scaleM;

    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 0);

    for (size_t i = 0; i < n_elem; i++) {
        const size_t *local_nodes = &elem_nodes[nl * i];
        SET_ELEM_INFO(nl, local_nodes, idx_map, coords, x_node, num);
        for (size_t q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            r = (model->m_type == AXISYMMETRIC) ? x_loc[0] : 1.;
            scaleM = w[q] * det * r;
            for (size_t j = 0; j < nl; j++) {
                size_t row = num[j];
                for (size_t k = 0; k < nl; k++) {
                    size_t col = num[k];
                    if (col <= row) { // only fill lower part
                        M->a[row][col] += phi[q][j] * phi[q][k] * scaleM;
                    }
                }
            }
        }
    }
    return;
}

void set_lsq_rhs(FE_Model *model, const double *u, double *rhs) {
    int n_elem = model->n_elem;
    int nl = model->n_local;
    const size_t *elem_nodes = model->elem_nodes;
    int n_node = model->n_node;
    const double *coords = model->coords;
    const size_t *idx_map = model->idx_map;

    size_t nq;
    size_t num[4] = {0};
    double w[4], xi[4][2], phi[4][4], dph[4][4][2], r, ir, h[4][4], sig[4];
    double x_node[4][2], u_node[8], x_loc[2], dphi[4][2], det, scale;

    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 0);
    set_hooke_matrix(1., model->nu, model->m_type, h);

    for (size_t i = 0; i < n_elem; i++) {
        const size_t *e_nodes = &elem_nodes[nl * i];
        SET_ELEM_INFO_U(nl, e_nodes, idx_map, coords, u, x_node, u_node, num);
        for (size_t q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            ir = (model->m_type == AXISYMMETRIC) ? det / x_loc[0] : 0.;
            compute_local_stress(nl, phi[q], dphi, u_node, h, ir, sig);
            r = (model->m_type == AXISYMMETRIC) ? x_loc[0] : 1.;
            scale = w[q] * r;
            for (size_t j = 0; j < nl; j++) {
                size_t row = num[j];
                rhs[0 * n_node + row] += sig[0] * phi[q][j] * scale; // s_11
                rhs[1 * n_node + row] += sig[1] * phi[q][j] * scale; // s_22
                rhs[2 * n_node + row] += sig[2] * phi[q][j] * scale; // s_12
                rhs[3 * n_node + row] += sig[3] * phi[q][j] * scale; // s_33
            }
        }
    }
    return;
}

/**
 * @brief Compute the nodal stress for the stress least squares problem
 * @param model Finite element model
 * @param u Solution displacement field (2 * n_node)
 * @param stresses Contains the nodal stresses on exit (n_node x 9)
 */
void compute_nodal_stress_lsq(FE_Model *mdl, const double *u, double *stress) {
    SymBandMatrix *M = mdl->M_scalar;
    size_t n_node = mdl->n_node;
    const size_t *idx_map = mdl->idx_map;
    const double nu = mdl->nu;

    if (!M) {
        M = allocate_sym_band_matrix(mdl->n_node, mdl->node_band);
        assemble_mass_lsq(mdl, M);
        sym_band_LDL(M->data, M->n, M->k);
    } else {
        M = mdl->M_scalar;
    }

    double *rhs = calloc(4 * n_node, sizeof(double));

    // Compute the element stresses
    set_lsq_rhs(mdl, u, rhs);

    solve_sym_band(M->data, M->n, M->k, rhs + 0 * n_node);
    for (size_t i = 0; i < n_node; i++)
        stress[i * 9 + 0] = rhs[0 * n_node + idx_map[i]];
    solve_sym_band(M->data, M->n, M->k, rhs + 1 * n_node);
    for (size_t i = 0; i < n_node; i++)
        stress[i * 9 + 4] = rhs[1 * n_node + idx_map[i]];
    solve_sym_band(M->data, M->n, M->k, rhs + 2 * n_node);
    for (size_t i = 0; i < n_node; i++)
        stress[i * 9 + 1] = stress[i * 9 + 3] = rhs[2 * n_node + idx_map[i]];

    if (mdl->m_type == AXISYMMETRIC) {
        solve_sym_band(M->data, M->n, M->k, rhs + 3 * n_node);
        for (size_t i = 0; i < n_node; i++)
            stress[i * 9 + 8] = rhs[3 * n_node + idx_map[i]];
    } else if (mdl->m_type == PLANE_STRAIN) {
        for (size_t i = 0; i < n_node; i++)
            stress[i * 9 + 8] = nu * (stress[i * 9 + 0] + stress[i * 9 + 4]);
    } else if (mdl->m_type == PLANE_STRESS) {
    } else {
        printf("Unknown model type: %d\n", mdl->m_type);
        exit(EXIT_FAILURE);
    }

    cblas_dscal(9 * n_node, mdl->E, stress, 1);
    free(rhs);
}

// s_node = 1/area int_patch s_elem dA
//        = 1/area sum_elem int_elem s_elem dA
//        = 1/area sum_elem sum_q w_q s_elem(xi_q) det(J)
void compute_nodal_stress_avg(FE_Model *mdl, const double *u, double *stress) {
    size_t nq, num[4];
    double w[4], xi[4][2], phi[4][4], dph[4][4][2];
    double det, dphi[4][2], x_loc[2], sig[4];
    double u_node[8], x_node[4][2], h[4][4], scale, ir;

    int n_elem = mdl->n_elem;
    int nl = mdl->n_local;
    const size_t *elem_nodes = mdl->elem_nodes;
    const double *coords = mdl->coords;
    const size_t *idx_map = mdl->idx_map;
    double nu = mdl->nu;

    set_hooke_matrix(1., nu, mdl->m_type, h);
    compute_shape_functions(mdl->e_type, &nq, w, xi, phi, dph, 0);
    double *den = calloc(mdl->n_node, sizeof(double));

    for (size_t i = 0; i < n_elem; i++) {
        const size_t *e_nodes = &elem_nodes[nl * i];
        SET_ELEM_INFO_U(nl, e_nodes, idx_map, coords, u, x_node, u_node, num);
        for (int q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            ir = (mdl->m_type == AXISYMMETRIC) ? det / x_loc[0] : 0.;
            compute_local_stress(nl, phi[q], dphi, u_node, h, ir, sig);
            for (size_t j = 0; j < nl; j++) {
                node_idx = e_nodes[j] - 1;
                den[node_idx] += w[q] * phi[q][j] * det;
                scale = w[q] * phi[q][j];
                // 1/det in eps/sig and integral det cancel out
                stress[9 * node_idx + 0] += scale * sig[0];
                stress[9 * node_idx + 1] += scale * sig[2];
                stress[9 * node_idx + 3] += scale * sig[2];
                stress[9 * node_idx + 4] += scale * sig[1];
                stress[9 * node_idx + 8] += scale * sig[3];
            }
        }
    }
    for (size_t j = 0; j < mdl->n_node; j++) {
        stress[9 * j + 0] *= mdl->E / den[j];
        stress[9 * j + 1] *= mdl->E / den[j];
        stress[9 * j + 3] *= mdl->E / den[j];
        stress[9 * j + 4] *= mdl->E / den[j];
        stress[9 * j + 8] *= mdl->E / den[j];
    }
    if (mdl->m_type == PLANE_STRAIN) {
        for (size_t j = 0; j < mdl->n_node; j++) {
            stress[9 * j + 8] = nu * (stress[9 * j + 0] + stress[9 * j + 4]);
        }
    }
    free(den);
}

/**
 * @brief Compute the forces on the boundary edges
 * @param model Finite element model
 * @param stress Stress tensor (n_node x (3*3))
 * @param data Array that stores (x, y, z, fx, fy, fz) for each edge
 * @param n_steps Number of steps
 * @param s Step (eigenmode)
 */
void compute_edge_forces(
    FE_Model *model, const double *stress, double *data, int n_steps, int s
) {
    size_t n1, n2;
    size_t idx_x, idx_f;
    size_t nnb = model->n_bd_edge;
    double norm, nx, ny, s11, s12, s22;

    for (size_t i = 0; i < nnb; i++) {
        idx_x = i * (3 + 3 * n_steps);
        idx_f = i * (3 + 3 * n_steps) + 3 + 3 * s;
        n1 = model->bd_edges[4 * i + 0] - 1;
        n2 = model->bd_edges[4 * i + 1] - 1;
        nx = model->coords[2 * n1 + 0] + model->coords[2 * n2 + 0];
        ny = model->coords[2 * n1 + 1] + model->coords[2 * n2 + 1];
        data[idx_x + 0] = nx / 2. * model->L_ref;
        data[idx_x + 1] = ny / 2. * model->L_ref;
        if (model->m_type == AXISYMMETRIC &&
            fabs(data[idx_x + 0]) < ZERO_RADIUS) {
            data[idx_x + 0] = 0.;
            continue;
        }
        nx = model->coords[2 * n2 + 1] - model->coords[2 * n1 + 1];
        ny = model->coords[2 * n1 + 0] - model->coords[2 * n2 + 0];
        norm = hypot(nx, ny);
        nx /= norm;
        ny /= norm;
        s11 = (stress[9 * n1 + 0] + stress[9 * n2 + 0]) / 2.;
        s12 = (stress[9 * n1 + 1] + stress[9 * n2 + 1]) / 2.;
        s22 = (stress[9 * n1 + 4] + stress[9 * n2 + 4]) / 2.;
        data[idx_f + 0] = s11 * nx + s12 * ny;
        data[idx_f + 1] = s12 * nx + s22 * ny;
        data[idx_f + 2] = 0.;
    }
}

// void compute_bd_forces(
//     FE_Model *model, const double *sol, double *data, int n_steps, int s
// ) {
//     size_t node_idx, nq, num, n1, n2, e, l;
//     double w[4], xi[4][2], phi[4][4], dph[4][4][2];
//     double det, dphi[4][2], x_loc[2], sig[4];
//     double u_node[8], x_node[4][2], h[4][4], norm, r, ir;

//     const size_t nl = model->n_local;
//     const size_t nnb = model->n_bd_edge;
//     const size_t *elem_nodes = model->elem_nodes;
//     const double *coords = model->coords;
//     const size_t *idx_map = model->idx_map;
//     int idx_x, idx_f;

//     compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 1);
//     double ref_l = cblas_dsum(nq, w, 1);
//     set_hooke_matrix(1., model->nu, model->m_type, h);

//     // With P1/Q1, the average value over the edge is the midpoint value
//     for (size_t i = 0; i < nnb; i++) {
//         idx_x = i * (3 + 3 * n_steps);
//         idx_f = i * (3 + 3 * n_steps) + 3 + 3 * s;
//         n1 = model->bd_edges[4 * i + 0] - 1;
//         n2 = model->bd_edges[4 * i + 1] - 1;
//         e = model->bd_edges[4 * i + 2];
//         l = model->bd_edges[4 * i + 3];
//         double nx = model->coords[2 * n2 + 1] - model->coords[2 * n1 + 1];
//         double ny = model->coords[2 * n1 + 0] - model->coords[2 * n2 + 0];
//         norm = hypot(nx, ny);
//         nx /= norm;
//         ny /= norm;
//         for (size_t j = 0; j < nl; j++) {
//             node_idx = elem_nodes[nl * e + (l + j) % nl] - 1;
//             num = idx_map[node_idx];
//             x_node[j][0] = coords[2 * node_idx + 0];
//             x_node[j][1] = coords[2 * node_idx + 1];
//             u_node[2 * j + 0] = sol[2 * num + 0];
//             u_node[2 * j + 1] = sol[2 * num + 1];
//         }
//         data[idx_x + 0] = (x_node[0][0] + x_node[1][0]) / 2. * model->L_ref;
//         data[idx_x + 1] = (x_node[0][1] + x_node[1][1]) / 2. * model->L_ref;
//         data[idx_x + 2] = 0.;
//         data[idx_f + 0] = data[idx_f + 1] = data[idx_f + 2] = 0.;
//         if (data[idx_x + 0] < ZERO_RADIUS && model->m_type == AXISYMMETRIC)
//             continue;
//         for (int q = 0; q < nq; q++) {
//             set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
//             r = (model->m_type == AXISYMMETRIC) ? x_loc[0] : 1.;
//             ir = (model->m_type == AXISYMMETRIC) ? 1. / x_loc[0] : 0.;
//             cblas_dscal(8, 1. / det, &dphi[0][0], 1);
//             compute_local_stress(nl, phi[q], dphi, u_node, h, ir, sig);
//             norm = w[q] * model->E * r / ref_l;
//             data[idx_f + 0] += norm * (sig[0] * nx + sig[2] * ny);
//             data[idx_f + 1] += norm * (sig[2] * nx + sig[1] * ny);
//             data[idx_f + 2] += 0.;
//         }
//     }
// }

double von_mises(double t[9]) {
    double v1 = SQUARE(t[0] - t[4]) + SQUARE(t[4] - t[8]) + SQUARE(t[8] - t[0]);
    double v2 = SQUARE(t[1]) + SQUARE(t[2]) + SQUARE(t[5]);
    return sqrt(0.5 * v1 + 3. * v2);
}

void tensor_eigws(double t[9], Model2D m_type, double eigws[2]) {
    // s_11  s_12  0
    // s_12  s_22  0
    //  0     0   s_33
    double tmp1, tmp2, l1, l2, l3;
    tmp1 = (t[0] + t[4]) / 2.;
    tmp2 = t[0] * t[4] - t[1] * t[1];
    l1 = tmp1 + sqrt(tmp1 * tmp1 - tmp2);
    l2 = tmp1 - sqrt(tmp1 * tmp1 - tmp2);
    l3 = t[8];
    eigws[0] = fmin(fmin(l1, l2), l3);
    eigws[1] = fmax(fmax(l1, l2), l3);
}

int compare_double(const void *a, const void *b) {
    return (*(double *)a > *(double *)b) - (*(double *)a < *(double *)b);
}

double compute_percentile(size_t n_node, const double *u_sol, double prc) {
    int idx;
    double norm;
    double *norms = malloc(n_node * sizeof(double));
    for (size_t i = 0; i < n_node; i++) {
        norm = hypot(u_sol[2 * i + 0], u_sol[2 * i + 1]);
        norms[i] = norm;
    }
    qsort(norms, n_node, sizeof(double), compare_double);
    idx = MAX(0, MIN(n_node - 1, (int)(prc * n_node)));
    return norms[idx];
}

/**
 * @brief Add the displacement field in Gmsh
 * @param model Finite element model
 * @param sol Solution vector (ux, uy)
 * @param views Gmsh views
 * @param n_steps Number of steps
 * @param s Step (eigenmodes)
 * @param force Array that stores the forces on the edges
 * @param bounds Bounds of the views
 */
void visualize_stress(
    FE_Model *model,
    const double *sol,
    const int *views,
    int n_steps,
    int s,
    double *force,
    double *bds
) {
    int ierr;
    size_t nn = model->n_node;
    size_t *nodes = malloc(nn * sizeof(size_t));
    double eig_stress[2];
    const int inds[4] = {0, 1, 4, 8};
    const int size = 9 * nn;
    double *sig_lsq = calloc(size, sizeof(double));
    compute_nodal_stress_lsq(model, sol, sig_lsq);
    compute_edge_forces(model, sig_lsq, force, n_steps, s);
    // compute_bd_forces(model, sol, force, n_steps, s);
#if DISPLAY_AVG
    double *sig_avg = calloc(size, sizeof(double));
    compute_nodal_stress_avg(model, sol, sig_avg);
#endif
    if (model->m_type != AXISYMMETRIC && CART_TO_POLAR) {
        cartesian_to_polar(nn, sig_lsq, model->coords);
    }
    for (size_t i = 0; i < nn; i++) {
        nodes[i] = i + 1;
        UPDATE_MAX(bds[2 * 2 + 1], von_mises(&sig_lsq[9 * i]));
        tensor_eigws(&sig_lsq[9 * i], model->m_type, eig_stress);
        UPDATE_MIN(bds[2 * 3 + 0], eig_stress[1]);
        UPDATE_MAX(bds[2 * 3 + 1], eig_stress[1]);
        UPDATE_MIN(bds[2 * 4 + 0], eig_stress[0]);
        UPDATE_MAX(bds[2 * 4 + 1], eig_stress[0]);
        for (int j = 0; j < 4; j++) {
            UPDATE_MIN(bds[2 * (5 + j) + 0], -fabs(sig_lsq[9 * i + inds[j]]));
            UPDATE_MAX(bds[2 * (5 + j) + 1], +fabs(sig_lsq[9 * i + inds[j]]));
        }
#if DISPLAY_AVG
        UPDATE_MAX(bds[2 * 9 + 1], von_mises(&sig_avg[9 * i]));
        tensor_eigws(&sig_avg[9 * i], model->m_type, eig_stress);
        UPDATE_MIN(bds[2 * 10 + 0], eig_stress[1]);
        UPDATE_MAX(bds[2 * 10 + 1], eig_stress[1]);
        UPDATE_MIN(bds[2 * 11 + 0], eig_stress[0]);
        UPDATE_MAX(bds[2 * 11 + 1], eig_stress[0]);
        for (int j = 0; j < 4; j++) {
            UPDATE_MIN(bds[2 * (12 + j) + 0], -fabs(sig_avg[9 * i + inds[j]]));
            UPDATE_MAX(bds[2 * (12 + j) + 1], +fabs(sig_avg[9 * i + inds[j]]));
        }
#endif
    }

    char *name;
    char *dtype = "NodeData";
    gmshModelGetCurrent(&name, &ierr);

    PLOT(views[2], s, name, dtype, nodes, nn, sig_lsq, size, 0., 9, -1, &ierr);
    free(sig_lsq);
#if DISPLAY_AVG
    PLOT(views[9], s, name, dtype, nodes, nn, sig_avg, size, 0., 9, -1, &ierr);
    free(sig_avg);
#endif
    free(nodes);
    free(name);
}

/**
 * @brief Add the displacement field in Gmsh
 * @param model Finite element model
 * @param sol Solution vector (ux, uy)
 * @param view Gmsh view_tag
 * @param step Step (eigenmodes)
 * @param bds Bounds of the views
 */
void visualize_disp(
    FE_Model *model, const double *sol, int view, int step, double *bounds
) {
    int ierr;
    size_t num;
    size_t nn = model->n_node;
    size_t *nodes = malloc(nn * sizeof(size_t));
    double *disp = malloc(3 * nn * sizeof(double));
    double norm;

    for (size_t i = 0; i < nn; i++) {
        nodes[i] = i + 1;
        num = model->idx_map[i];
        disp[3 * i + 0] = sol[2 * num + 0] * model->L_ref;
        disp[3 * i + 1] = sol[2 * num + 1] * model->L_ref;
        disp[3 * i + 2] = 0.;
        norm = hypot(disp[3 * i + 0], disp[3 * i + 1]);
        bounds[1] = fmax(bounds[1], norm);
    }
    bounds[0] = fmax(bounds[0], compute_percentile(nn, disp, 0.95));

    char *name;
    char *dtype = "NodeData";
    gmshModelGetCurrent(&name, &ierr);
    PLOT(view, step, name, dtype, nodes, nn, disp, 3 * nn, 0., 3, -1, &ierr);
    gmshViewOptionSetNumber(view, "VectorType", 5, &ierr);
    gmshViewOptionSetNumber(view, "DrawPoints", 0, &ierr);
    free(nodes);
    free(disp);
    free(name);
}

/**
 * @brief Add the boundary forces in Gmsh
 * @param model Finite element model
 * @param data List of size n_views * n_bd_edge * 6
 * @param view List of views
 * @param n_steps Number of steps in gmsh
 * @param bounds List of bounds associated to each view
 * @note Needs all steps to be computed (by visualize_stress)
 */
void visualize_bd_forces(
    FE_Model *model, const double *data, int view, int n_steps, double *bds
) {
    int ierr;
    char *name;
    double norm, max_force;
    size_t nnb = model->n_bd_edge;
    size_t idx;
    int incx = 3 + 3 * n_steps;

    for (size_t s = 0; s < n_steps; s++) {
        max_force = 0.;
        for (size_t i = 0; i < nnb; i++) {
            idx = i * (incx) + (3 + 3 * s);
            norm = hypot(data[idx + 0], data[idx + 1]);
            max_force = fmax(max_force, norm);
        }
        bds[1] = fmax(bds[1], max_force);
    }

    gmshModelGetCurrent(&name, &ierr);
    gmshViewAddListData(view, "VP", nnb, data, nnb * incx, &ierr);
    gmshViewOptionSetNumber(view, "VectorType", 2, &ierr);
    free(name);
}

/**
 * @brief Add the Gmsh views
 * @param views_ptr Pointer to the Gmsh view list
 * @param n_views_ptr Pointer to the number of views
 * @param bounds_ptr Pointer to the bounds of the views
 */
void add_gmsh_views(int **views_ptr, int *n_views_ptr, double **bounds_ptr) {
    int ierr, *prev_views;
    size_t n_prev;
    gmshViewGetTags(&prev_views, &n_prev, &ierr);
    for (int i = 0; i < n_prev; i++) {
        gmshViewRemove(prev_views[i], &ierr);
    }
    free(prev_views);
    int n_views = 2 + 7 + 7 * DISPLAY_AVG;
    int *views = malloc(n_views * sizeof(int));
    double *bounds = malloc(2 * n_views * sizeof(double));
    for (int i = 0; i < n_views; i++) {
        bounds[2 * i + 0] = bounds[2 * i + 1] = 0.;
        views[i] = -1;
    }
    views[0] = gmshViewAdd("forces", -1, &ierr);       // Boundary forces
    views[1] = gmshViewAdd("displacement", -1, &ierr); // Displacement
    views[2] = gmshViewAdd("stress lsq", -1, &ierr);   // Stress Least Squares
#if DISPLAY_AVG
    views[9] = gmshViewAdd("stress avg", -1, &ierr); // Stress Average
#endif

    *views_ptr = views;
    *n_views_ptr = n_views;
    *bounds_ptr = bounds;
    return;
}

void create_tensor_aliases(int *views) {
    int ierr;
    views[3] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress max eigw
    gmshViewOptionSetNumber(views[3], "TensorType", 2, &ierr);
    views[4] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress max eigw
    gmshViewOptionSetNumber(views[4], "TensorType", 3, &ierr);

    views[5] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress xx
    gmshViewOptionSetNumber(views[5], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[5], "ComponentMap0", 0, &ierr);
    views[6] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress xy
    gmshViewOptionSetNumber(views[6], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[6], "ComponentMap0", 1, &ierr);
    views[7] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress yy
    gmshViewOptionSetNumber(views[7], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[7], "ComponentMap0", 4, &ierr);
    views[8] = gmshViewAddAlias(views[2], 1, -1, &ierr); // Stress zz
    gmshViewOptionSetNumber(views[8], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[8], "ComponentMap0", 8, &ierr);

#if DISPLAY_AVG
    views[10] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress max eigw
    gmshViewOptionSetNumber(views[10], "TensorType", 2, &ierr);
    views[11] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress min eigw
    gmshViewOptionSetNumber(views[11], "TensorType", 3, &ierr);

    views[12] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress xx
    gmshViewOptionSetNumber(views[12], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[12], "ComponentMap0", 0, &ierr);
    views[13] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress xy
    gmshViewOptionSetNumber(views[13], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[13], "ComponentMap0", 1, &ierr);
    views[14] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress yy
    gmshViewOptionSetNumber(views[14], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[14], "ComponentMap0", 4, &ierr);
    views[15] = gmshViewAddAlias(views[9], 1, -1, &ierr); // Stress zz
    gmshViewOptionSetNumber(views[15], "ForceNumComponents", 1, &ierr);
    gmshViewOptionSetNumber(views[15], "ComponentMap0", 8, &ierr);
#endif
}

/**
 * @brief Set the bounds of the views
 * @param n_views Number of views
 * @param views Gmsh view list
 * @param bounds Bounds of the views (min1, max1, min2, max2, ...)
 * @param mode 1: Linear system solve, 2: Eigenmodes
 */
void set_view_options(int n_views, int *views, double *bounds) {
    int ierr;

    // Set view field ranges
    double min, max;
    for (int i = 0; i < n_views; i++) {
        min = bounds[2 * i + 0];
        max = bounds[2 * i + 1];
        gmshViewOptionSetNumber(views[i], "RangeType", 2, &ierr);
        gmshViewOptionSetNumber(views[i], "CustomMin", min, &ierr);
        gmshViewOptionSetNumber(views[i], "CustomMax", max, &ierr);
        // gmshViewOptionSetNumber(views[i], "ShowScale", 0, &ierr);
    }
    gmshViewOptionSetNumber(views[1], "RangeType", 2, &ierr);
    gmshViewOptionSetNumber(views[1], "CustomMin", 0., &ierr);

    // Set displacement factor
    double factor, dl, box[6];
    gmshModelGetBoundingBox(-1, -1, UNZIP6(box), &ierr);
    dl = hypot(box[3] - box[0], box[4] - box[1]);
    // factor = bounds[2 * 1 + 1];
    factor = bounds[2 * 1 + 0];
    factor = (factor < 1e-20) ? 1. : (dl / 20.) / factor;
    gmshViewOptionSetNumber(views[1], "DisplacementFactor", factor, &ierr);

    // Hide the mesh
    gmshOptionSetNumber("Mesh.SurfaceEdges", 0, &ierr);

    // Hide the fields
    for (int i = 0; i < n_views; i++)
        gmshViewOptionSetNumber(views[i], "Visible", 0, &ierr);
    gmshViewOptionSetNumber(views[1], "Visible", 1, &ierr);
}

/**
 * @brief Revolve the geometry in Gmsh
 * @param model Finite element model
 */
void revolve_geometry(FE_Model *model) {
    if (model->m_type != AXISYMMETRIC)
        return;
    int ierr;
    int *dim_tags, *out_dim_tags;
    size_t dim_tags_n, out_dim_tags_n;
    gmshModelGetEntities(&dim_tags, &dim_tags_n, 2, &ierr);
    printf("Dim tags: %zu\n", dim_tags_n);
    for (size_t i = 0; i < dim_tags_n; i++) {
        printf("Dim tag: %d\n", dim_tags[i]);
    }

    // gmshModelOccCopy(
    //     dim_tags, dim_tags_n, &out_dim_tags, &out_dim_tags_n, &ierr
    // );
    // gmshModelOccRotate(
    //     out_dim_tags, out_dim_tags_n, 0., 0., 0., 0., 1., 0., 1e-1, &ierr
    // );

    // clang-format off
    gmshModelOccRevolve(
        dim_tags, dim_tags_n, 0., 0., 0., 0., 1., 0., 2.0 * M_PI, 
        &out_dim_tags, &out_dim_tags_n, NULL, 0, NULL, 0, 0, &ierr
    );
    // clang-format on

    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Geometry.Surfaces", 1, &ierr);
    gmshOptionSetNumber("Geometry.SurfaceType", 1, &ierr);
    free(dim_tags);
    free(out_dim_tags);
}