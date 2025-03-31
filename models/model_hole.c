#include "models.h"
#include "utils_models.h"
#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// Rectangular plate with a hole
// -----------------------------------------------------------------------------

static double scale;
static const double W = 2.0;
static const double H = 2.0;
static const double R = 0.1;

static const double _E = 211e9;
static const double _nu = 0.30;
static const double _rho = 7.85e3;
static const double _L = R;

void set_physics_hole(double params[4], Model2D *type) {
    params[0] = _E;
    params[1] = _nu;
    params[2] = _rho;
    params[3] = _L;
    *type = PLANE_STRESS;
}

void set_bk_source_hole(double rho, const double xy[2], double f[2]) {
    double x = xy[0] * _L;
    double y = xy[1] * _L;
    f[0] = 0. * x * y;
    f[1] = 0.;
}

void set_bd_disp_hole(int entity, char kind, const double xy[2], double u[1]) {
    u[0] = 0.;
}

void set_bd_force_hole(int entity, char kind, const double xy[2], double f[2]) {
    if (entity == 4 && kind == 'x') {
        f[0] = 0.;
        f[1] = 1e4;
    } else {
        f[0] = 0.;
        f[1] = 0.;
    }
}

double size_field_hole(
    int dim, int tag, double x, double y, double z, double lc, void *data
) {
    double h1 = 0.02;
    double h2 = 0.20;
    double dref = 0.5;
    double d = hypot(x, y) - R;
    return scale * hermite(d, dref, h1, h2);
}

void mesh_hole(double mesh_size_factor, int e_type) {
    int ierr;
    int idRect =
        gmshModelOccAddRectangle(0., 0., 0., W / 2., H / 2., -1, 0.0, &ierr);
    int idHole =
        gmshModelOccAddDisk(0., 0., 0., R, R, -1, NULL, 0, NULL, 0, &ierr);
    int objs[2] = {2, idRect};
    int tool[2] = {2, idHole};
    gmshModelOccCut(
        objs, 2, tool, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr
    );
    gmshModelOccSynchronize(&ierr);

    int force_x[] = {4};
    int force_y[] = {};
    int force_n[] = {};
    int force_t[] = {};
    int dirichlet_x[] = {2};
    int dirichlet_y[] = {5};
    int dirichlet_n[] = {};
    int dirichlet_t[] = {};
    gmshModelAddPhysicalGroup(1, force_x, 1, 1, "force_x", &ierr);
    gmshModelAddPhysicalGroup(1, force_y, 0, 2, "force_y", &ierr);
    gmshModelAddPhysicalGroup(1, force_n, 0, 3, "force_n", &ierr);
    gmshModelAddPhysicalGroup(1, force_t, 0, 4, "force_t", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_x, 1, 5, "fix_x", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_y, 1, 6, "fix_y", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_n, 0, 7, "fix_n", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_t, 0, 8, "fix_t", &ierr);

    scale = mesh_size_factor; // set global variable
    gmshModelMeshSetSizeCallback(size_field_hole, NULL, &ierr);
    ierr = mesh_tri_quad(e_type);
    // gmshFltkRun(&ierr);

    return;
}
