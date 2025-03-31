#include "models.h"
#include "utils_models.h"
#include <FL/math.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// Annulus (cylinder slice) under uniform inner / outer pressure
// -----------------------------------------------------------------------------

static double scale;
static const double Ri = 0.1;
static const double Ro = 0.15;
static const double Pi = 5e3;
static const double Po = 1e3;

static const double _E = 40e9;
static const double _nu = 0.20;
static const double _rho = 2300.;
static const double _L = Ri;

void set_physics_tank_cut(double params[4], Model2D *type) {
    params[0] = _E;
    params[1] = _nu;
    params[2] = _rho;
    params[3] = _L;
    *type = PLANE_STRAIN;
}

void set_bk_source_tank_cut(double rho, const double xy[2], double f[2]) {
    double x = xy[0] * _L;
    double y = xy[1] * _L;
    f[0] = 0. * x * y;
    f[1] = 0.;
}

void set_bd_disp_tank_cut(
    int entity, char kind, const double xy[2], double u[1]
) {
    u[0] = 0.;
}

void set_bd_force_tank_cut(
    int entity, char kind, const double xy[2], double f[2]
) {
    if (entity == 5 && kind == 'n') {
        f[0] = 0.;
        f[1] = -Po;
    } else if (entity == 7 && kind == 'n') {
        f[0] = 0.;
        f[1] = +Pi;
    } else {
        f[0] = 0.;
        f[1] = 0.;
    }
}

double size_field_tank_cut(
    int dim, int tag, double x, double y, double z, double lc, void *data
) {
    double d, h1, h2;
    double dr = Ro - Ri;
    double hmin = dr / 4.;
    double hmax = dr / 3.;
    double dref = dr / 3.;
    d = Ro - hypot(x, y);
    h1 = scale * hermite(d, dref, hmin, hmax);
    d = hypot(x, y) - Ri;
    h2 = scale * hermite(d, dref, hmin, hmax);
    return fmin(h1, h2);
}

void mesh_tank_cut(double mesh_size_factor, int e_type) {

    compute_stress_tank(Ri, Ro, Pi, Po, _nu);

    int ierr;
    int p1 = gmshModelOccAddPoint(Ri, 0., 0., 0., -1, &ierr);
    int p2 = gmshModelOccAddPoint(Ro, 0., 0., 0., -1, &ierr);
    int p3 = gmshModelOccAddPoint(0., Ro, 0., 0., -1, &ierr);
    int p4 = gmshModelOccAddPoint(0., Ri, 0., 0., -1, &ierr);
    int l1 = gmshModelOccAddLine(p1, p2, -1, &ierr);
    int l2 = gmshModelOccAddCircle(
        0., 0., 0., Ro, -1, 0., M_PI / 2., NULL, 0, NULL, 0, &ierr
    );
    int l3 = gmshModelOccAddLine(p3, p4, -1, &ierr);
    int l4 = gmshModelOccAddCircle(
        0., 0., 0., Ri, -1, 0., M_PI / 2., NULL, 0, NULL, 0, &ierr
    );
    int curves[4] = {l1, l2, l3, l4};
    gmshModelOccAddCurveLoop(curves, 4, -1, &ierr);
    int s1 = gmshModelOccAddPlaneSurface(&l1, 1, -1, &ierr);
    int remove[] = {1, l1, 1, l2, 1, l3, 1, l4};
    gmshModelOccRemove(remove, 8, 0, &ierr);
    gmshModelOccSynchronize(&ierr);

    int force_x[] = {};
    int force_y[] = {};
    int force_n[] = {5, 7};
    int force_t[] = {};
    int dirichlet_x[] = {6};
    int dirichlet_y[] = {1};
    int dirichlet_n[] = {};
    int dirichlet_t[] = {};
    gmshModelAddPhysicalGroup(1, force_x, 0, 1, "force_x", &ierr);
    gmshModelAddPhysicalGroup(1, force_y, 0, 2, "force_y", &ierr);
    gmshModelAddPhysicalGroup(1, force_n, 2, 3, "force_n", &ierr);
    gmshModelAddPhysicalGroup(1, force_t, 0, 4, "force_t", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_x, 1, 5, "fix_x", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_y, 1, 6, "fix_y", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_n, 0, 7, "fix_n", &ierr);
    gmshModelAddPhysicalGroup(1, dirichlet_t, 0, 8, "fix_t", &ierr);
    gmshModelAddPhysicalGroup(2, &s1, 1, 9, "tank", &ierr);

    scale = mesh_size_factor; // set global variable
    gmshModelMeshSetSizeCallback(size_field_tank_cut, NULL, &ierr);
    gmshModelMeshSetSizeFromBoundary(2, s1, 0, &ierr);
    ierr = mesh_tri_quad(e_type);

    // gmshFltkRun(&ierr);
    // exit(0);

    return;
}
