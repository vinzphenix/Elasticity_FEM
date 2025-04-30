/**
 * File:    model_tank_axi.c
 * Author:  Vincent Degrooff
 * Created: 2025
 * Description: 
 *   Set the physics / geometry / boundary conditions of 
 *   a thick-walled cylinder under pressure (axisymmetric formulation)
 * Project: FEM Simulation Toolkit for Linear Elasticity
 */

#include "models.h"
#include "utils_models.h"
#include <math.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double scale;
static const double Ri = 0.1;
static const double Ro = 0.15;
static const double Pi = 5e3;
static const double Po = 1e3;

static const double _E = 40e9;
static const double _nu = 0.20;
static const double _rho = 2300.;
static const double _L = Ro;
static const double H = 0.2;

void set_physics_tank_axi(double params[4], Model2D *type) {
    params[0] = _E;
    params[1] = _nu;
    params[2] = _rho;
    params[3] = _L;
    *type = AXISYMMETRIC;
}

void set_bk_source_tank_axi(double rho, const double xy[2], double f[2]) {
    double r = xy[0] * _L;
    double z = xy[1] * _L;
    f[0] = 0. + 0. * r * z;
    f[1] = 0.;
}

void set_bd_disp_tank_axi(
    int entity, char kind, const double xy[2], double u[1]
) {
    u[0] = 0.;
}

void set_bd_force_tank_axi(
    int entity, char kind, const double xy[2], double f[2]
) {
    if (entity == 2 && kind == 'n') {
        f[0] = 0.;
        f[1] = -Po;
    } else if (entity == 4 && kind == 'n') {
        f[0] = 0.;
        f[1] = -Pi;
    } else {
        f[0] = 0.;
        f[1] = 0.;
    }
}

double size_field_tank_axi(
    int dim, int tag, double x, double y, double z, double lc, void *data
) {
    double d, h1, h2;
    double dr = Ro - Ri;
    double hmin = dr / 4.;
    double hmax = dr / 3.;
    double dref = dr / 3.;
    d = Ro - x;
    h1 = scale * hermite(d, dref, hmin, hmax);
    d = x - Ri;
    h2 = scale * hermite(d, dref, hmin, hmax);
    return fmin(h1, h2);
}

void mesh_tank_axi(double mesh_size_factor, int e_type) {

    // compute_stress_tank(Ri, Ro, Pi, Po, _nu);

    int ierr;
    int rect = gmshModelOccAddRectangle(Ri, 0., 0., Ro - Ri, H, -1, 0.0, &ierr);
    gmshModelOccSynchronize(&ierr);
    int force_x[] = {};
    int force_y[] = {};
    int force_n[] = {2, 4};
    int force_t[] = {};
    int dirichlet_x[] = {};
    int dirichlet_y[] = {1, 3};
    int dirichlet_n[] = {};
    int dirichlet_t[] = {};
    gmshModelAddPhysicalGroup(1, force_x, 0, 1, "force_x", &ierr);
    gmshModelAddPhysicalGroup(1, force_y, 0, 2, "force_y", &ierr);
    gmshModelAddPhysicalGroup(1, force_n, 2, 3, "force_n", &ierr);
    gmshModelAddPhysicalGroup(1, force_t, 0, 4, "force_t", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_x, 0, 5, "fix_x", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_y, 2, 6, "fix_y", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_n, 0, 7, "fix_n", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_t, 0, 8, "fix_t", &ierr);
    scale = mesh_size_factor; // set global variable

    // gmshModelMeshSetTransfiniteSurface(1, "AlternateLeft", NULL, 0, &ierr);

    gmshModelMeshSetSizeCallback(size_field_tank_axi, NULL, &ierr);
    gmshModelMeshSetSizeFromBoundary(2, rect, 0, &ierr);
    ierr = mesh_tri_quad(e_type);
    // gmshFltkRun(&ierr);
    // exit(0);

    return;
}
