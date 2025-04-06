#include "model.h"
#include "matrix.h"
#include "models.h"
#include "renumber.h"

#include <gmshc.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

FE_Model *create_FE_Model(
    const char *name, ElementType etp, Renumbering rnb, LinearSolver solver
) {

    FE_Model *model = (FE_Model *)malloc(sizeof(FE_Model));
    model->model_name = name;

    model->e_type = etp;    // element: triangle, quadrilateral
    model->renum = rnb;     // renumbering: none, x, y, rcmk
    model->solver = solver; // solver: Band, CG_[flavor]
    model->M_scalar = NULL;
    double parameters[4];

    if (strcmp(name, "beam") == 0) {
        model->mesh_model = mesh_beam;
        model->set_bk_source = set_bk_source_beam;
        model->set_bd_force = set_bd_force_beam;
        model->set_bd_disp = set_bd_disp_beam;
        set_physics_beam(parameters, &model->m_type);
    } else if (strcmp(name, "exam") == 0) {
        model->mesh_model = mesh_exam;
        model->set_bk_source = set_bk_source_exam;
        model->set_bd_force = set_bd_force_exam;
        model->set_bd_disp = set_bd_disp_exam;
        set_physics_exam(parameters, &model->m_type);
    } else if (strcmp(name, "hole") == 0) {
        model->mesh_model = mesh_hole;
        model->set_bk_source = set_bk_source_hole;
        model->set_bd_force = set_bd_force_hole;
        model->set_bd_disp = set_bd_disp_hole;
        set_physics_hole(parameters, &model->m_type);
    } else if (strcmp(name, "section") == 0) {
        model->mesh_model = mesh_section;
        model->set_bk_source = set_bk_source_section;
        model->set_bd_force = set_bd_force_section;
        model->set_bd_disp = set_bd_disp_section;
        set_physics_section(parameters, &model->m_type);
    } else if (strcmp(name, "fork") == 0) {
        model->mesh_model = mesh_fork;
        model->set_bk_source = set_bk_source_fork;
        model->set_bd_force = set_bd_force_fork;
        model->set_bd_disp = set_bd_disp_fork;
        set_physics_fork(parameters, &model->m_type);
    } else if (strcmp(name, "tank_axi") == 0) {
        model->mesh_model = mesh_tank_axi;
        model->set_bk_source = set_bk_source_tank_axi;
        model->set_bd_force = set_bd_force_tank_axi;
        model->set_bd_disp = set_bd_disp_tank_axi;
        set_physics_tank_axi(parameters, &model->m_type);
    } else if (strcmp(name, "tank_cut") == 0) {
        model->mesh_model = mesh_tank_cut;
        model->set_bk_source = set_bk_source_tank_cut;
        model->set_bd_force = set_bd_force_tank_cut;
        model->set_bd_disp = set_bd_disp_tank_cut;
        set_physics_tank_cut(parameters, &model->m_type);
    } else {
        printf("Unknown model: %s\n", name);
        exit(EXIT_FAILURE);
    }

    model->E = parameters[0];     // Young's modulus
    model->nu = parameters[1];    // Poisson's ratio
    model->rho = parameters[2];   // Density
    model->L_ref = parameters[3]; // Reference length (just for scaling)

    return model;
}

void free_FE_Model(FE_Model *model) {
    free(model->elem_nodes);
    free(model->coords);
    free(model->idx_map);
    free(model->bd_edges);
    free_band_sym(model->M);
    free_band_sym(model->K);
    if (model->M_scalar != NULL) {
        free_band_sym(model->M_scalar);
    }
    free(model);
}

void load_mesh(FE_Model *model) {
    int ierr;
    int e_type = model->e_type;
    gmshModelMeshRebuildNodeCache(1, &ierr);
    double *coord_bad, *coords;
    size_t *elem_tags, *node_tags, *e_nodes;
    size_t tmp, max_node, n_node, n_elem;

    gmshModelMeshGetNodes(
        &node_tags, &n_node, &coord_bad, &tmp, NULL, NULL, 2, -1, 1, 0, &ierr
    );
    gmshModelMeshGetElementsByType(
        e_type, &elem_tags, &n_elem, &e_nodes, &tmp, -1, 0, 1, &ierr
    );
    max_node = 0;
    for (size_t i = 0; i < n_node; i++) {
        max_node = MAX(max_node, node_tags[i]);
    }
    if (max_node != n_node) {
        printf("Nodes not contig. : max: %zu  nn: %zu\n", max_node, n_node);
        exit(EXIT_FAILURE);
    }
    coords = (double *)malloc(2 * n_node * sizeof(double));
    for (size_t i = 0; i < n_node; i++) {
        tmp = node_tags[i] - 1; // 0-based indexing
        coords[2 * tmp + 0] = coord_bad[3 * i + 0] / model->L_ref;
        coords[2 * tmp + 1] = coord_bad[3 * i + 1] / model->L_ref;
    }

    model->n_elem = n_elem;
    model->n_local = e_type + 1;
    model->elem_nodes = e_nodes;
    model->n_node = n_node;
    model->coords = coords;
    model->e_tags = elem_tags;

    find_boundary_nodes(model);

    free(coord_bad);
    free(node_tags);
}

/**
 * @brief Renumber the nodes of a mesh
 * @param model Finite element model
 * @param renum_type Renumbering type (1: x, 2: y, 3: RCMK)
 */
void renumber_nodes(FE_Model *model) {
    renumber(
        model->n_elem,
        model->n_local,
        model->elem_nodes,
        model->n_bd_edge,
        model->bd_edges,
        model->n_node,
        model->coords,
        &model->idx_map,
        model->renum
    );
    model->node_band = compute_band(
        model->n_elem, model->n_local, model->elem_nodes, model->idx_map
    );
}

int sort_edge(const void *a, const void *b) {
    Edge *edge_a = (Edge *)a;
    Edge *edge_b = (Edge *)b;
    size_t a_min = MIN(edge_a->n1, edge_a->n2);
    size_t b_min = MIN(edge_b->n1, edge_b->n2);
    if (a_min != b_min) {
        return a_min - b_min;
    } else {
        size_t a_max = MAX(edge_a->n1, edge_a->n2);
        size_t b_max = MAX(edge_b->n1, edge_b->n2);
        return a_max - b_max;
    }
}

void edge_compress(Edge *edges, size_t n_edges, size_t *ne, size_t *n_bd_e) {
    size_t n = 0;
    size_t n_bd = 0;
    for (size_t i = 0; i < n_edges; i++, n++) {
        size_t min1 = MIN(edges[i + 0].n1, edges[i + 0].n2);
        size_t max1 = MAX(edges[i + 0].n1, edges[i + 0].n2);
        size_t min2 = MIN(edges[i + 1].n1, edges[i + 1].n2);
        size_t max2 = MAX(edges[i + 1].n1, edges[i + 1].n2);
        // No problem at i = n_edges - 1 because of ghost edge
        edges[n] = edges[i];
        if (i < n_edges - 1 && min1 == min2 && max1 == max2) {
            edges[n].e2 = edges[i + 1].e1;
            i++;
        } else { // Boundary edge
            n_bd++;
        }
    }
    *ne = n;
    *n_bd_e = n_bd;
}

/**
 * @brief Compute the edge list of a mesh
 * @param n_loc Number of nodes per element
 * @param n_elem Number of elements
 * @param elems Element nodes (1-based indexing)
 * @param n_bd_edge_ptr Number of edges on the boundary
 * @param n_edge_ptr Number of unique edges
 * @param edges_ptr List of unique edges
 */
void create_edges(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t *n_bd_edge_ptr,
    size_t *n_edge_ptr,
    Edge **edges_ptr
) {
    size_t s = 0; // keep 1-based indexing
    Edge *edges = malloc((n_loc * n_elem + 1) * sizeof(Edge)); // 1 ghost edge
    for (size_t i = 0; i < n_elem; i++) {
        for (size_t j = 0; j < n_loc; j++) {
            edges[n_loc * i + j].e1 = i;
            edges[n_loc * i + j].e2 = n_elem; // Mark as boundary edge
            edges[n_loc * i + j].n1 = elems[n_loc * i + j] - s;
            edges[n_loc * i + j].n2 = elems[n_loc * i + (j + 1) % n_loc] - s;
        }
    }

    // Sort edges
    qsort(edges, n_loc * n_elem, sizeof(Edge), sort_edge);

    // Compress edges
    size_t n_edges, n_bd;
    edge_compress(edges, n_loc * n_elem, &n_edges, &n_bd);
    edges = realloc(edges, n_edges * sizeof(Edge));

    *n_bd_edge_ptr = n_bd;
    *n_edge_ptr = n_edges;
    *edges_ptr = edges;
}

int get_local_edge(size_t n1, size_t n2, size_t n_loc, size_t *nodes) {
    for (size_t i = 0; i < n_loc; i++) {
        if (nodes[i] == n1 && nodes[(i + 1) % n_loc] == n2) {
            return i;
        }
    }
    return -1;
}

/**
 * @brief Find all boundary nodes in a mesh
 * @param model Finite element model
 * @note The boundary nodes are stored in model->bd_nodes
 */
void find_boundary_nodes(FE_Model *model) {

    size_t n_bd_edge, n_edge;
    Edge *edges;
    size_t n_loc = model->n_local;
    size_t n_elem = model->n_elem;
    size_t *elems = model->elem_nodes;
    create_edges(n_loc, n_elem, elems, &n_bd_edge, &n_edge, &edges);

    // Allocate boundary edges
    const int n_info = 4;
    model->n_bd_edge = n_bd_edge;
    model->bd_edges = malloc(n_info * n_bd_edge * sizeof(size_t));
    for (size_t i = 0, j = 0; i < n_edge; i++) {
        if (edges[i].e2 == n_elem) {
            model->bd_edges[n_info * j + 0] = edges[i].n1;
            model->bd_edges[n_info * j + 1] = edges[i].n2;
            model->bd_edges[n_info * j + 2] = edges[i].e1;
            model->bd_edges[n_info * j + 3] = get_local_edge(
                edges[i].n1, edges[i].n2, n_loc, &elems[n_loc * edges[i].e1]
            );
            j++;
        }
    }

    free(edges);
}
