#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "model.h"
#include <stddef.h>

#define SET_ELEM_INFO(nl, e_nodes, idx_map, coords, x_node, num)               \
    size_t node_idx;                                                           \
    for (size_t j = 0; j < (nl); j++) {                                        \
        node_idx = (e_nodes)[j] - 1;                                           \
        (x_node)[j][0] = (coords)[2 * node_idx + 0];                           \
        (x_node)[j][1] = (coords)[2 * node_idx + 1];                           \
        (num)[j] = (idx_map)[node_idx];                                        \
    }

#define SET_ELEM_INFO_U(nl, e_nodes, idx_map, coords, u, x_node, u_node, num)  \
    size_t node_idx;                                                           \
    for (size_t j = 0; j < (nl); j++) {                                        \
        node_idx = (e_nodes)[j] - 1;                                           \
        (x_node)[j][0] = (coords)[2 * node_idx + 0];                           \
        (x_node)[j][1] = (coords)[2 * node_idx + 1];                           \
        (num)[j] = (idx_map)[node_idx];                                        \
        (u_node)[2 * j + 0] = (u)[2 * (num)[j] + 0];                           \
        (u_node)[2 * j + 1] = (u)[2 * (num)[j] + 1];                           \
    }

void compute_shape_functions(
    ElementType e_type,
    size_t *nq,
    double w[4],
    double xi[4][2],
    double phi[4][4],
    double dph[4][4][2],
    const int boundary
);
void set_geometry(
    const int n_loc,
    const double x_node[4][2],
    const double phi[4],
    const double dph[4][2],
    double x_ptr[2],
    double det[1],
    double dphi[4][2]
);
void set_strain_basis(
    const int n_loc,
    const double phi[4],
    const double dph[4][2],
    const double inv_r,
    double B[8][4]
);
void set_hooke_matrix(double E, double nu, int m_type, double h[4][4]);

void assemble_system(FE_Model *model);
void add_bulk_source(FE_Model *model, double *rhs);
void enforce_bd_conditions(FE_Model *model, double *rhs);

#endif
