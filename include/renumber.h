#ifndef RENUMBER_H
#define RENUMBER_H
#include <stddef.h>

typedef struct {
    size_t n_node;
    size_t max_degree;
    size_t *offsets;   // offset[i] points to the adjacency of node i
    size_t *adjacency; // flattened array of neighbors
} MeshGraph;

typedef struct {
    size_t node;
    size_t degree;
} NodeDegree;

MeshGraph *create_adjacency(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t n_node,
    size_t n_bd_edge,
    size_t *bd_edges
);

void free_mesh_graph(MeshGraph *graph);

void rcmk_renumber(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t n_node,
    size_t n_bd_node,
    size_t *bd_nodes,
    size_t *idx_map
);

void renumber(
    size_t n_elem,
    size_t n_loc,
    size_t *elems,
    size_t n_bd_node,
    size_t *bd_nodes,
    size_t n_node,
    double *coords,
    size_t **idx_map,
    int kind
);

size_t compute_band(
    const size_t ne, const size_t nl, const size_t *elem, const size_t *idx_map
);

#endif