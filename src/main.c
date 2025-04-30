#define _POSIX_C_SOURCE 199309L
#include "elasticity.h"
#include "model.h"
#include "power.h"
#include "visualize.h"
#include <FL/math.h>
#include <cblas.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERBOSE 1
#define PRECISION 10
#define TIMESP(t1, t2)                                                         \
    ((double)((t2).tv_sec - (t1).tv_sec) +                                     \
     1e-9 * ((double)((t2).tv_nsec - (t1).tv_nsec)))

void save_solutions(
    FE_Model *model, char *name, double *eigv, double *eigw, int nb
) {
    FILE *f = fopen(name, "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening file %s\n", name);
        exit(1);
    }
    size_t nn = model->n_node;
    int idx;

    // Write solution(s)
    nb = MAX(nb, 1);
    fprintf(f, "# %d\n", nb);
    for (int i_eig = 0; i_eig < nb; i_eig++) {
        fprintf(f, "# %zu ", nn);
        if (eigw != NULL) 
            fprintf(f, "%20.15le", eigw[i_eig]);
        fprintf(f, "\n");
        for (int i = 0; i < nn; i++) {
            idx = i_eig * 2 * nn + 2 * model->idx_map[i];
            fprintf(f, "%.*le ", PRECISION, eigv[idx + 0]);
            fprintf(f, "%.*le\n", PRECISION, eigv[idx + 1]);
        }
    }
    fclose(f);
}

void solve_deformation(FE_Model *model) {
    struct timespec t1, t2, t3;
    size_t n_node = model->n_node;
    SymBandMatrix *K = model->K;
    double *rhs = (double *)calloc(2 * n_node, sizeof(*rhs));
    double *sol = (double *)malloc(2 * n_node * sizeof(*sol));

    clock_gettime(CLOCK_MONOTONIC, &t1);
    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);
    // write_band_sym(model->K, rhs, "K.txt");
    
    clock_gettime(CLOCK_MONOTONIC, &t2);
    printf("%30s : %6.3lf s\n", "Boundaries", TIMESP(t1, t2));

    int nit = solve_system(K, model->solver, rhs, sol);
    clock_gettime(CLOCK_MONOTONIC, &t3);
    printf("%30s : %6.3lf s", "Sys solve", TIMESP(t2, t3));
    if (1 != nit)
        printf(" (%d it.)\n", nit);
    printf("\n\n");

    int ierr, n_views, *views;
    double *bounds;
    add_gmsh_views(&views, &n_views, &bounds);

    double *data_forces = malloc(6 * model->n_bd_edge * sizeof(double));
    // compute_bd_forces(model, rhs, data_forces, 1, 0);
    visualize_disp(model, sol, views[1], 0, &bounds[2]);
    visualize_stress(model, sol, views, 1, 0, data_forces, bounds);
    visualize_bd_forces(model, data_forces, views[0], 1, &bounds[0]);

    create_tensor_aliases(views);
    set_view_options(n_views, views, bounds);
    // revolve_geometry(model);
    gmshFltkRun(&ierr);
    gmshFltkFinalize(&ierr);
    free(data_forces);
    free(rhs);
    free(sol);
    free(bounds);
    free(views);
}

void find_eigenmodes(FE_Model *model, int nb) {
    int ierr;
    struct timespec t1, t2;
    int n = model->K->n;
    int max_it = 10000;
    int n_bd_edge = model->n_bd_edge;
    double *eigw = (double *)malloc(nb * (n + 1) * sizeof(double));
    double *eigv = eigw + nb;
    double rtol = 1e-11;

    double *rhs = (double *)malloc(2 * n * sizeof(double));
    enforce_bd_conditions(model, rhs);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    compute_eigvs_deflation(model->K, model->M, eigw, eigv, nb, rtol, max_it);
    clock_gettime(CLOCK_MONOTONIC, &t2);
    printf("\n%17s%2d %s : %6.3lf s\n", "", nb, "Eigenmodes", TIMESP(t1, t2));
    for (int i = 0; i < nb; i++) {
        eigw[i] = fmax(eigw[i], 0.);
        eigw[i] = sqrt(eigw[i]) * sqrt(model->E / model->rho) / model->L_ref;
        eigw[i] /= (2. * M_PI);
    }
    // save_solutions(model, "eigenmodes.msh", eigv, eigw, nb);
    // print_XAX(model->K->data, &eigv[0], n, model->K->k, nb, "X' K X");
    // print_XAX(model->M->data, &eigv[0], n, model->K->k, nb, "X' M X");
    // print_XAX(NULL, &eigv[0], n, model->K->k, nb, "X' X");

    int n_views, *views;
    double *bounds;
    add_gmsh_views(&views, &n_views, &bounds);
    double *forces = malloc(n_bd_edge * (3 + 3 * nb) * sizeof(double));

    printf("\n");
    for (int i = 0; i < nb; i++) {
        visualize_disp(model, &eigv[i * n], views[1], i, &bounds[2]);
        visualize_stress(model, &eigv[i * n], views, nb, i, forces, bounds);
        printf("%26s %3d : %8.2lf Hz\n", "Eigenmode", i + 1, eigw[i]);
    }
    printf("\n");

    visualize_bd_forces(model, forces, views[0], nb, bounds);
    create_tensor_aliases(views);
    set_view_options(n_views, views, bounds);
    // revolve_geometry(model);
    gmshFltkRun(&ierr);
    gmshFltkFinalize(&ierr);
    free(forces);
    free(eigw);
    free(rhs);
    free(bounds);
    free(views);
}

void display_info(FE_Model *model, int step, struct timespec ts[4]) {

    char *m_str[3] = {"Plane stress", "Plane strain", "Axisymmetric"};
    char *r_str[4] = {"No", "X", "Y", "RCMK"};

    if (step == 1) {
        printf(
            "\n===========  Linear elasticity simulation - FEM  ===========\n\n"
        );
        printf("%30s = %s\n", "Model", model->model_name);
        printf("%30s = %s\n", "Model type", m_str[model->m_type]);
        printf("%30s = %.3e\n", "Young's Modulus E", model->E);
        printf("%30s = %.3e\n", "Poisson ratio nu", model->nu);
        printf("%30s = %.3e\n\n", "Density rho", model->rho);
    } else if (step == 2) {
        char *e_str = (model->e_type == TRI) ? "Triangle" : "Quadrilateral";
        printf("%30s = %s\n", "Element type", e_str);
        printf("%30s = %zu\n", "Number of elements", model->n_elem);
        printf("%30s = %zu\n", "Number of nodes", model->n_node);
        printf("%30s = %s\n", "Renumbering", r_str[model->renum]);
        printf("%30s = %zu\n", "Matrix bandwidth", 2 * model->node_band + 1);
        const char *s_str = solver_name(model->solver);
        printf("%30s = %s\n\n", "Solver", s_str);
    } else if (step == 3) {
        printf("%30s : %6.3lf s\n", "Mesh load", TIMESP(ts[0], ts[1]));
        printf("%30s : %6.3lf s\n", "Node renumber", TIMESP(ts[1], ts[2]));
        printf("%30s : %6.3lf s\n", "Assembly", TIMESP(ts[2], ts[3]));
    }
}

int main(int argc, char *argv[]) {

    int ierr, nb_mode;
    double meshSizeFactor;
    struct timespec times[4];
    if ((argc < 4) || (sscanf(argv[2], "%d", &nb_mode)) != 1 ||
        (sscanf(argv[3], "%lf", &meshSizeFactor)) != 1) {
        printf("Usage: \n./deformation <model> <nb modes> <meshSizeFactor>\n");
        printf("model: one of the model implemented in models/\n");
        printf("nb modes: number of eigenmodes OR 0 for static loading\n");
        printf("meshSizeFactor: mesh size factor for gmsh\n");
        return -1;
    }

    // Simulation parameters
    const ElementType e_type = QUAD;
    const Renumbering renum = RENUM_RCMK;
    const LinearSolver solver = CG_ILU0;

    FE_Model *model = create_FE_Model(argv[1], e_type, renum, solver);
    display_info(model, 1, NULL);

    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 2, &ierr);
    model->mesh_model(meshSizeFactor, e_type);

    clock_gettime(CLOCK_MONOTONIC, &times[0]);
    load_mesh(model);

    clock_gettime(CLOCK_MONOTONIC, &times[1]);
    renumber_nodes(model);
    display_info(model, 2, NULL);

    clock_gettime(CLOCK_MONOTONIC, &times[2]);
    assemble_system(model);
    clock_gettime(CLOCK_MONOTONIC, &times[3]);
    display_info(model, 3, times);

    if (nb_mode == 0)
        solve_deformation(model);
    else
        find_eigenmodes(model, nb_mode);

    // Free stuff
    gmshFinalize(&ierr);
    free_FE_Model(model);
    return 0;
}
