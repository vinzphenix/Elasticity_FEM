#ifndef VISUALIZE_H
#define VISUALIZE_H

#include "model.h"
#include <stddef.h>

void visualize_disp(
    FE_Model *model, const double *sol, int view, int step, double *bounds
);
void visualize_stress(
    FE_Model *model,
    const double *sol,
    const int *views,
    int n_steps,
    int s,
    double *force,
    double *bds
);
void visualize_bd_forces(
    FE_Model *model, const double *data, int view, int n_steps, double *bds
);

void compute_nodal_stress_lsq(FE_Model *mdl, const double *u, double *stress);
void compute_nodal_stress_avg(FE_Model *mdl, const double *u, double *stress);
void cartesian_to_polar(
    size_t n_node, const size_t *num, const double *x, double *s, double *u
);

void add_gmsh_views(int **views_ptr, int *n_views_ptr, double **bounds_ptr);
void set_view_options(int n_views, int *views, double *bounds);
void create_tensor_aliases(int *views);
void revolve_geometry(FE_Model *model);

#endif