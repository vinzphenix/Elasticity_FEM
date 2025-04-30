/**
 * File:    renumber.c
 * Author:  Vincent Degrooff
 * Created: 2025
 *
 * Description:
 *   Node renumbering strategies
 * 
 * Project:
 *   FEM Simulation Toolkit for Linear Elasticity
 */

#include "renumber.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

void print_adjacency(MeshGraph *graph) {
    printf("Number of nodes   = %3zu\n", graph->n_node);
    printf("Maximum adjacency = %3zu\n", graph->max_degree);
    for (size_t i = 0; i < graph->n_node; i++) {
        printf("Adjacency of node %3zu : ", i + 1);
        for (size_t j = graph->offsets[i]; j < graph->offsets[i + 1]; j++) {
            printf(" %3zu, ", graph->adjacency[j] + 1);
        }
        printf("\n");
    }
}

void add_neighbor(size_t *offsets, size_t *adjacency, size_t n1, size_t n2) {
    size_t ptr = offsets[n1];
    adjacency[ptr] = n2;
    offsets[n1] += 1;
}

/**
 * @brief Create the adjacency of a mesh
 * @param n_loc Number of nodes per element
 * @param n_elem Number of elements
 * @param elems Element nodes (1-based indexing)
 * @param n_node Number of nodes
 * @param n_bd_edge Number of boundary edges
 * @param bd_edges Boundary edges (1-based indexing)
 * @return MeshGraph Adjacency structure
 */
MeshGraph *create_adjacency(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t n_node,
    size_t n_bd_edge,
    size_t *bd_edges
) {

    // Create adjacency graph
    const size_t s = 1; // 1-based indexing
    MeshGraph *graph = malloc(sizeof(MeshGraph));
    graph->n_node = n_node;
    graph->offsets = calloc(n_node + 1, sizeof(size_t));

    // First pass: count degree of each node
    // Shift the count by 1 in the "offset" array
    size_t node;
    for (size_t i = 0; i < n_elem; i++) {
        for (size_t j = 0; j < n_loc; j++) {
            node = elems[n_loc * i + j] - s;
            graph->offsets[node + 1] += 1;
        }
    }
    // Add 1 for every boundary node valence
    for (size_t i = 0; i < n_bd_edge; i++) {
        node = bd_edges[4 * i + 0] - s;
        graph->offsets[node + 1] += 1;
    }

    // Compute required size of adjacency
    graph->max_degree = 0;
    for (size_t i = 1; i <= n_node; i++) {
        graph->max_degree = MAX(graph->max_degree, graph->offsets[i]);
        graph->offsets[i] = graph->offsets[i] + graph->offsets[i - 1];
    }
    graph->adjacency = calloc(graph->offsets[n_node], sizeof(size_t));

    // Second pass: add mapping
    size_t n1, n2;
    for (size_t i = 0; i < n_elem; i++) {
        for (size_t j = 0; j < n_loc; j++) { // loop over edges
            n1 = elems[i * n_loc + j] - s;
            n2 = elems[i * n_loc + (j + 1) % n_loc] - s;
            add_neighbor(graph->offsets, graph->adjacency, n1, n2);
        }
    }
    // Add missing connections on boundary
    for (size_t i = 0; i < n_bd_edge; i++) {
        n1 = bd_edges[4 * i + 0] - s;
        n2 = bd_edges[4 * i + 1] - s;
        add_neighbor(graph->offsets, graph->adjacency, n2, n1);
    }

    // Shift back the offset array
    for (size_t i = n_node; i > 0; i--) {
        graph->offsets[i] = graph->offsets[i - 1];
    }
    graph->offsets[0] = 0;

    // print_adjacency(graph);
    return graph;
}

void free_mesh_graph(MeshGraph *graph) {
    free(graph->adjacency);
    free(graph->offsets);
    free(graph);
}

int compare_by_degree(const void *a, const void *b) {
    return ((NodeDegree *)a)->degree - ((NodeDegree *)b)->degree;
}

/**
 * @brief Reverse Cuthill-McKee renumbering
 * @param n_loc Number of nodes per element
 * @param n_elem Number of elements
 * @param elems Element nodes (1-based indexing)
 * @param n_node Number of nodes
 * @param idx_map Mapping between original and optimized indices
 */
void rcmk_renumber(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t n_node,
    size_t n_bd_edge,
    size_t *bd_edges,
    size_t *idx_map
) {
    MeshGraph *graph =
        create_adjacency(n_loc, n_elem, elems, n_node, n_bd_edge, bd_edges);

    size_t *offsets = graph->offsets;
    size_t *adjacency = graph->adjacency;

    size_t *bfs_nodes = malloc(n_node * sizeof(size_t));
    size_t *visited = idx_map;
    memset(visited, 0, n_node * sizeof(size_t));
    NodeDegree *neighs = malloc(graph->max_degree * sizeof(NodeDegree));

    size_t start, node, neigh, n_neigh, front, rear;
    start = 0;
    n_neigh = n_node - 1;
    // Pick peripheral node with the smallest degree
    for (size_t i = 0; i < n_bd_edge; i++) {
        node = bd_edges[4 * i + 0] - 1;
        size_t degree = offsets[node + 1] - offsets[node];
        if (degree < n_neigh) {
            n_neigh = degree;
            start = node;
        }
    }

    front = 0, rear = 0;
    bfs_nodes[rear++] = start;
    visited[start] = 1;

    while (front < rear) {
        node = bfs_nodes[front++];
        n_neigh = 0;
        for (size_t j = offsets[node]; j < offsets[node + 1]; j++) {
            neigh = adjacency[j];
            if (!visited[neigh]) {
                neighs[n_neigh].degree = offsets[neigh + 1] - offsets[neigh];
                neighs[n_neigh].node = neigh; // add to current adjacency
                visited[neigh] = 1;           // mark as visited
                n_neigh++;
            }
        }
        qsort(neighs, n_neigh, sizeof(NodeDegree), compare_by_degree);
        for (size_t j = 0; j < n_neigh; j++) {
            bfs_nodes[rear++] = neighs[j].node; // put the node in the set
        }
    }

    // Renumber (with reversed)
    for (size_t i = 0; i < n_node; i++) {
        idx_map[bfs_nodes[i]] = n_node - 1 - i;
    }

    free(neighs);
    free(bfs_nodes);
    free_mesh_graph(graph);
}

/**
 * @brief Compute the band size of a linear system from a mesh
 * @param ne Number of elements
 * @param nl Number of nodes per element
 * @param elem Element nodes (1-based indexing)
 * @param idx_map Mapping between original and optimized indices
 * @return size_t Maximum neighborhood distance
 */
size_t compute_band(
    const size_t ne, const size_t nl, const size_t *elem, const size_t *idx_map
) {
    size_t idx, imax, imin;
    size_t band = 0;
    size_t s = 1; // 1-based indexing
    idx = 0;
    for (size_t e = 0; e < ne; e++) {
        imin = imax = idx_map[elem[idx + 0] - s];
        for (size_t i = 1; i < nl; i++) {
            imin = MIN(imin, idx_map[elem[idx + i] - s]);
            imax = MAX(imax, idx_map[elem[idx + i] - s]);
        }
        band = MAX(band, imax - imin);
        idx += nl;
    }
    return band;
}

double *cmp_XY;
int compare_node_x(const void *n1, const void *n2) {
    return (cmp_XY[2 * (*(size_t *)n1) + 0] > cmp_XY[2 * (*(size_t *)n2) + 0]) -
           (cmp_XY[2 * (*(size_t *)n1) + 0] < cmp_XY[2 * (*(size_t *)n2) + 0]);
}
int compare_node_y(const void *n1, const void *n2) {
    return (cmp_XY[2 * (*(size_t *)n1) + 1] > cmp_XY[2 * (*(size_t *)n2) + 1]) -
           (cmp_XY[2 * (*(size_t *)n1) + 1] < cmp_XY[2 * (*(size_t *)n2) + 1]);
}

/**
 * @brief Renumber nodes in a mesh
 * @param n_elem Number of elements
 * @param n_loc Number of nodes per element
 * @param elems Element nodes (1-based indexing)
 * @param n_node Number of nodes
 * @param coords Node coordinates (x, y)
 * @param idx_map Mapping between original and optimized indices
 * @param kind Renumbering kind: 1: x, 2: y, 3: RCMK
 */
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
) {
    *idx_map = malloc(n_node * sizeof(size_t));
    if (kind == 0) {
        for (size_t i = 0; i < n_node; i++) {
            (*idx_map)[i] = i;
        }
    } else if (kind == 1 || kind == 2) {
        cmp_XY = coords;
        size_t *indices = malloc(n_node * sizeof(size_t));
        for (size_t i = 0; i < n_node; i++) {
            indices[i] = i;
        }
        if (kind == 1) {
            qsort(indices, n_node, sizeof(size_t), compare_node_x);
        } else {
            qsort(indices, n_node, sizeof(size_t), compare_node_y);
        }
        for (size_t i = 0; i < n_node; i++) {
            (*idx_map)[indices[i]] = i;
        }
        free(indices);
    } else if (kind == 3) {
        rcmk_renumber(
            n_loc, n_elem, elems, n_node, n_bd_node, bd_nodes, *idx_map
        );
    } else {
        fprintf(stderr, "Unknown renumbering kind %d\n", kind);
        exit(1);
    }
}
