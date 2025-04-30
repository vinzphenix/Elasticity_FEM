/**
 * File:    conv_beam.c
 * Author:  Vincent Degrooff
 * Created: 2025
 * 
 * Description: 
 *   Runs a convergence study on the beam problem.
 * 
 * Project: 
 *   FEM Simulation Toolkit for Linear Elasticity
 */

#include "elasticity.h"
#include "model.h"
#include "power.h"
#include "visualize.h"
#include <cblas.h>
#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SQ(x) ((x) * (x))
#define M_PI 3.14159265358979323846

static const double H = 1.0;  // HARDCODED
static const double G = 9.81; // HARDCODED
static const double gam[4] = {
    1.87510406871, 4.69409113297, 7.85475743824, 10.9955407349
};

double get_beam_solution(double x, int mode, FE_Model *model) {
    double L = model->L_ref;
    double E = model->E;
    double rho = model->rho;
    double xi = x / L;

    double res;
    if (mode == 0) {
        res = -rho * G * L * L * L * L / (2. * E * H * H);
        res *= xi * xi * (6. - 4. * xi + xi * xi);
    } else {
        res = gam[mode - 1];
        res = (sin(res) - sinh(res)) * (sin(res * xi) - sinh(res * xi)) +
              (cos(res) + cosh(res)) * (cos(res * xi) - cosh(res * xi));
    }
    return res;
}

void get_neutral_axis_sol(
    FE_Model *model,
    int mode,
    size_t n_pts,
    const double *sol,
    double *cut,
    double *ana
) {
    int ierr, n_views, *views;
    double *bounds;
    add_gmsh_views(&views, &n_views, &bounds);
    visualize_disp(model, sol, views[1], 0, &bounds[2]);

    double factor = (mode == 0) ? 10. : 0.01;
    gmshViewOptionSetNumber(views[1], "DisplacementFactor", factor, &ierr);

    gmshPluginSetString("CutParametric", "X", "u", &ierr);
    gmshPluginSetString("CutParametric", "Y", "0", &ierr);
    gmshPluginSetString("CutParametric", "Z", "0", &ierr);
    gmshPluginSetNumber("CutParametric", "MinU", 0., &ierr);
    gmshPluginSetNumber("CutParametric", "MaxU", model->L_ref, &ierr);
    gmshPluginSetNumber("CutParametric", "NumPointsU", n_pts, &ierr);
    gmshPluginSetNumber("CutParametric", "NumPointsV", 1, &ierr);
    gmshPluginSetNumber("CutParametric", "View", 1, &ierr);
    gmshPluginSetNumber("CutParametric", "ConnectPoints", 0, &ierr);
    int cut_view = gmshPluginRun("CutParametric", &ierr);

    char **data_types;
    size_t data_types_n, num_elem_n, *data_n, data_nn;
    int *num_elem;
    double **data;
    // clang-format off
    gmshViewGetListData(
        cut_view, &data_types, &data_types_n, &num_elem, &num_elem_n, 
        &data, &data_n, &data_nn, 0, &ierr
    );
    // clang-format on

    double x1, v1, v1e, norm, sign;
    for (int i = 0; i < num_elem[0]; i++) {
        x1 = data[0][6 * i + 0];
        v1 = data[0][6 * i + 4];
        v1e = get_beam_solution(x1, mode, model);
        cut[i] = v1;
        ana[i] = v1e;
    }

    if (0 < mode) {
        sign = (0. < cut[n_pts - 1]) ? 1. : -1.;
        norm = cblas_dnrm2(n_pts, cut, 1);
        cblas_dscal(n_pts, sign / norm, cut, 1);
        sign = (0. < ana[n_pts - 1]) ? 1. : -1.;
        norm = cblas_dnrm2(n_pts, ana, 1);
        cblas_dscal(n_pts, sign / norm, ana, 1);
    }
    
    gmshFltkRun(&ierr);
    free(bounds);
    free(views);
    // should free gmsh allocated memory...
}

void deformation(
    FE_Model *model, size_t n_pts, double *ref, int save_fine, double *errs
) {
    size_t n_node = model->n_node;
    SymBandMatrix *K = model->K;
    double *rhs = (double *)calloc(2 * n_node, sizeof(*rhs));
    double *sol = (double *)malloc(2 * n_node * sizeof(*sol));

    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);
    solve_system(K, model->solver, rhs, sol);

    double norm_ana, norm_ref, end_ana, end_ref;
    double *ana = (double *)malloc(n_pts * sizeof(double));
    double *cut = (double *)malloc(n_pts * sizeof(double));
    
    get_neutral_axis_sol(model, 0, n_pts, sol, cut, ana);
    if (save_fine)
        memcpy(ref, cut, n_pts * sizeof(double));

    end_ana = ana[n_pts - 1];
    end_ref = ref[n_pts - 1];
    norm_ana = cblas_dnrm2(n_pts, ana, 1);
    norm_ref = cblas_dnrm2(n_pts, ref, 1);
    cblas_daxpy(n_pts, -1., cut, 1, ana, 1);
    cblas_daxpy(n_pts, -1., ref, 1, cut, 1);
    errs[0] = cblas_dnrm2(n_pts, cut, 1) / norm_ref;
    errs[1] = fabs(cut[n_pts - 1] / end_ref);
    errs[2] = cblas_dnrm2(n_pts, ana, 1) / norm_ana;
    errs[3] = fabs(ana[n_pts - 1] / end_ana);

    free(cut);
    free(ana);
    free(rhs);
    free(sol);
    return;
}

void eigenmodes(
    FE_Model *model, int nb, size_t n_pts, double *ref, int save, double *errs
) {

    const double L = model->L_ref;
    const double E = model->E;
    const double rho = model->rho;
    size_t n = 2 * model->n_node;

    double **refs = (double **)malloc(nb * sizeof(double *));
    double **errors = (double **)malloc(nb * sizeof(double *));
    for (int i = 0; i < nb; i++) {
        refs[i] = ref + i * (1 + n_pts);
        errors[i] = errs + i * 4;
    }

    double *rhs = (double *)calloc(n, sizeof(*rhs));
    double *eigw = (double *)malloc(nb * (n + 1) * sizeof(double));
    double *eigv = eigw + nb;

    int max_it = 10000;
    double rtol = 1e-11;

    enforce_bd_conditions(model, rhs);
    compute_eigvs_deflation(model->K, model->M, eigw, eigv, nb, rtol, max_it);

    double om_ana, om_ref;
    double *ana = (double *)malloc(n_pts * sizeof(double));
    double *cut = (double *)malloc(n_pts * sizeof(double));

    for (int i = 0; i < nb; i++) {
        eigw[i] = fmax(eigw[i], 0.);
        eigw[i] = sqrt(eigw[i] * model->E / model->rho) / model->L_ref;
        
        get_neutral_axis_sol(model, i + 1, n_pts, eigv + i * n, cut, ana);
        if (save) {
            memcpy(refs[i], cut, n_pts * sizeof(double));
            refs[i][n_pts] = eigw[i];
        }

        om_ana = (SQ(gam[i] / L) * H * sqrt(E / 12. / rho));
        om_ref = refs[i][n_pts];

        printf("w[%d] = %9.4lf ? %9.4lf ? %9.4lf rad/s\n", i, eigw[i], om_ref, om_ana);

        cblas_daxpy(n_pts, -1., cut, 1, ana, 1);
        cblas_daxpy(n_pts, -1., refs[i], 1, cut, 1);
        errors[i][0] = fabs(eigw[i] - om_ref) / om_ref;
        errors[i][1] = cblas_dnrm2(n_pts, cut, 1);
        errors[i][2] = fabs(eigw[i] - om_ana) / om_ana;
        errors[i][3] = cblas_dnrm2(n_pts, ana, 1);
    }

    free(cut);
    free(ana);
    free(rhs);
    free(eigw);
    free(refs);
    free(errors);
}

// #define SAVE
void save_errors(FE_Model *model, int n_mode, double msf, double *errors) {
    char filename[64] = "";
    int inc = 4;

    strcat(filename, "./analysis/conv_beam");
    strcat(filename, n_mode == 0 ? "_static.txt" : "_eigen.txt");
    FILE *fp = fopen(filename, "a");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }
    // clang-format off
    if (n_mode == 0) {
#ifdef SAVE
        fprintf(fp, "%d %8.4lf %6zu %15.8le %15.8le %15.8le %15.8le\n",
            model->e_type, msf, model->n_node,
            errors[0], errors[1], errors[2], errors[3]
        );
#endif
        printf(
            "Mesh size factor: %8.4lf | %6zu nodes  | "
            "err = %9.3le, %9.3le, %9.3le, %9.3le\n",
            msf, model->n_node, errors[0], errors[1], errors[2], errors[3]
        );
    } else {
        for (int i = 0; i < n_mode; i++) {
#ifdef SAVE
            fprintf(
                fp, "%4d %8.4lf %6zu %4d %15.8le %15.8le %15.8le %15.8le\n", 
                model->e_type, msf, model->n_node, i + 1,
                errors[i * inc + 0], errors[i * inc + 1], 
                errors[i * inc + 2], errors[i * inc + 3]
            );
#endif
            printf(
                "Mesh size factor: %8.4lf | %6zu nodes  | mode %d | "
                "err = %9.3le, %9.3le, %9.3le, %9.3le\n",
                msf, model->n_node, i + 1,
                errors[i * inc + 0], errors[i * inc + 1],
                errors[i * inc + 2], errors[i * inc + 3]
            );
        }
    }
    fclose(fp);
    // clang-format on
}

int main(int argc, char *argv[]) {

    int ierr;
    double mesh_size_factor;
    FE_Model *model;
    size_t n_pts = 1000;

    // Simulation parameters
    const ElementType e_type = TRI;
    const Renumbering renum = RENUM_RCMK;
    const LinearSolver solver = Band;

    const int nb_mode = 3;
    double errors[4 * nb_mode];
    double *ref = (double *)malloc(nb_mode * (1 + n_pts) * sizeof(double));

    for (int s = 0; s < 6; s++) {
        mesh_size_factor = 0.20 * pow(2., s);
        model = create_FE_Model("beam", e_type, renum, solver);

        gmshInitialize(argc, argv, 0, 0, &ierr);
        gmshOptionSetNumber("General.Verbosity", 2, &ierr);
        model->mesh_model(mesh_size_factor, e_type);
        load_mesh(model);
        renumber_nodes(model);
        assemble_system(model);

        deformation(model, n_pts, ref, s == 0, errors);
        save_errors(model, 0, mesh_size_factor, errors);
        
        // eigenmodes(model, nb_mode, n_pts, ref, s == 0, errors);
        // save_errors(model, nb_mode, mesh_size_factor, errors);

        gmshFinalize(&ierr);
        free_FE_Model(model);
    }

    free(ref);
    return 0;
}
