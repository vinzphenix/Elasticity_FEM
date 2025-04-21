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

// ! HARDCODED
static const double Ri = 0.10;
static const double Ro = 0.15;
static const double Pi = 5e3;
static const double Po = 1e3;

double exact_solution(double E, double nu, double r, double sigma[3]) {
    double ri2 = Ri * Ri;
    double ro2 = Ro * Ro;
    double c1 = (Pi * ri2 - Po * ro2) / (ro2 - ri2);
    double c2 = ro2 * ri2 / (ro2 - ri2) * (Pi - Po);
    sigma[0] = c1 - c2 / (r * r);
    sigma[1] = c1 + c2 / (r * r);
    sigma[2] = 2. * nu * c1;
    return (1. + nu) / E * ((1. - 2. * nu) * c1 * r + c2 / r);
}

void compute_error(FE_Model *model, double *sol_u, double error[7]) {
    size_t nq;
    double w[4], xi[4][2], phi[4][4], dph[4][4][2];
    double x_node[4][2], dphi[4][2], x_loc[2], det;
    double r, ur_h, ur_e, sigma_e[3], sigma_lsq[3], sigma_avg[3], exact[7];

    compute_shape_functions(model->e_type, &nq, w, xi, phi, dph, 0);
    
    size_t nn = model->n_node;
    size_t num[4] = {0};
    size_t nl = model->n_local;
    int axisym = model->m_type == AXISYMMETRIC;
    double L_ref = model->L_ref;
    double *coords = model->coords;
    const size_t *idx_map = model->idx_map;
    error[0] = error[1] = error[2] = error[3] = 0.;
    exact[0] = exact[1] = exact[2] = exact[3] = 0.;
    
    int idx_s[3] = {0};
    idx_s[1] = (axisym) ? 8 : 4;
    idx_s[2] = (axisym) ? 4 : 8;
    
    double *stress_lsq = calloc(9 * nn, sizeof(double));
    double *stress_avg = calloc(9 * nn, sizeof(double));
    compute_nodal_stress_lsq(model, sol_u, stress_lsq);
    compute_nodal_stress_avg(model, sol_u, stress_avg);
    if (!axisym) {
        cartesian_to_polar(nn, idx_map, model->coords, stress_lsq, sol_u);
        cartesian_to_polar(nn, idx_map, model->coords, stress_avg, NULL);
    }

    for (size_t i = 0; i < model->n_elem; i++) {
        const size_t *local_nodes = &(model->elem_nodes[nl * i]);
        SET_ELEM_INFO(nl, local_nodes, idx_map, coords, x_node, num);
        for (size_t q = 0; q < nq; q++) {
            set_geometry(nl, x_node, phi[q], dph[q], x_loc, &det, dphi);
            r = (axisym) ? x_loc[0] : hypot(x_loc[0], x_loc[1]);
            ur_e = exact_solution(model->E, model->nu, r * L_ref, sigma_e);
            r = (axisym) ? x_loc[0] : 1.;
            ur_h = 0.;
            sigma_lsq[0] = sigma_lsq[1] = sigma_lsq[2] = 0.;
            sigma_avg[0] = sigma_avg[1] = sigma_avg[2] = 0.;
            for (size_t j = 0; j < nl; j++) {
                size_t jj = local_nodes[j] - 1;
                ur_h += sol_u[2 * num[j] + 0] * phi[q][j];
                for (size_t k = 0; k < 3; k++) {
                    sigma_lsq[k] += stress_lsq[9 * jj + idx_s[k]] * phi[q][j];
                    sigma_avg[k] += stress_avg[9 * jj + idx_s[k]] * phi[q][j];
                }
            }
            det = w[q] *r * det;
            exact[0] += det * SQUARE(ur_e);
            error[0] += det * SQUARE(ur_e - ur_h * L_ref);
            for (size_t k = 0; k < 3; k++) {
                exact[1+k] += det * SQUARE(sigma_e[k]);
                error[1+k] += det * SQUARE(sigma_e[k] - sigma_lsq[k]);
                error[4+k] += det * SQUARE(sigma_e[k] - sigma_avg[k]);
            }
            // printf("%15.5le  %15.5le | %15.5le  %15.5le | %15.5le  %15.5le\n", sigma_e[0], sigma_h[0], sigma_e[1], sigma_h[1], sigma_e[2], sigma_h[2]);
        }
    }

    error[0] = sqrt(error[0] / exact[0]);
    for (size_t k = 0; k < 3; k++) {
        error[1 + k] = sqrt(error[1 + k] / exact[1 + k]);
        error[4 + k] = sqrt(error[4 + k] / exact[1 + k]);
    }

    free(stress_lsq);
    free(stress_avg);
    return;
}

void solve_defo(FE_Model *model, double *sol) {
    size_t n_node = model->n_node;
    SymBandMatrix *K = model->K;
    double *rhs = (double *)calloc(2 * n_node, sizeof(*rhs));
    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);
    solve_system(K, model->solver, rhs, sol);
    free(rhs);
    return;
}

// #define SAVE
void save_errors(FE_Model *model, double msf, double *errors) {
    const char *m_type = (model->m_type == AXISYMMETRIC) ? "AXISYM" : "PLANE";
    char filename[64] = "./analysis/conv_tank.txt";
    FILE *fp = fopen(filename, "a");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }
    // clang-format off
#ifdef SAVE
    fprintf(
        fp, 
        "%5d %5d %8.4lf %6zu %15.8le "
        "%15.8le %15.8le %15.8le "
        "%15.8le %15.8le %15.8le\n",
        model->m_type, model->e_type, msf, model->n_node, errors[0], 
        errors[1], errors[2], errors[3], errors[4], errors[5], errors[6]
    );
#endif
    printf(
        "Mesh size factor: %8.4lf | %6zu nodes | %s | "
        "err = %9.3le, %9.3le, %9.3le, %9.3le, %9.3le, %9.3le, %9.3le\n",
        msf, model->n_node, m_type, errors[0], 
        errors[1], errors[2], errors[3], errors[4], errors[5], errors[6]
    );
    fclose(fp);
    // clang-format on
}

int main(int argc, char *argv[]) {

    int ierr;
    double mesh_size_factor, *sol;
    FE_Model *model;

    // Simulation parameters
    const ElementType e_type = QUAD;
    const Renumbering renum = RENUM_RCMK;
    const LinearSolver solver = Band;

    double errors[7];

    for (int s = 0; s < 5; s++) {
        mesh_size_factor = 0.025 * pow(2., s);
        // model = create_FE_Model("tank_cut", e_type, renum, solver);
        model = create_FE_Model("tank_axi", e_type, renum, solver);

        gmshInitialize(argc, argv, 0, 0, &ierr);
        gmshOptionSetNumber("General.Verbosity", 2, &ierr);
        model->mesh_model(mesh_size_factor, e_type);
        load_mesh(model);
        renumber_nodes(model);
        assemble_system(model);

        sol = (double *)malloc(2 * model->n_node * sizeof(*sol));
        solve_defo(model, sol);
        compute_error(model, sol, errors);
        save_errors(model, mesh_size_factor, errors);

        gmshFinalize(&ierr);
        free_FE_Model(model);
        free(sol);
    }

    return 0;
}
