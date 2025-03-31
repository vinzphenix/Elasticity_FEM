#ifndef MODEL_H
#define MODEL_H

#include "matrix.h"
#include "models.h"

#include <stddef.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define SQUARE(X) ((X) * (X))

typedef enum Renumbering {
    RENUM_NO,
    RENUM_X,
    RENUM_Y,
    RENUM_RCMK
} Renumbering;

typedef enum ElementType {
    TRI = 2,
    QUAD = 3
} ElementType;

typedef struct FE_Model {
    const char *model_name;
    double E;
    double nu;
    double rho;
    double L_ref;
    Model2D m_type;
    ElementType e_type;
    Renumbering renum;
    size_t node_band;
    size_t n_elem;
    size_t n_local;
    size_t n_node;
    size_t n_bd_edge;
    size_t *e_tags;     // 1-based indexing
    size_t *elem_nodes; // 1-based indexing
    size_t *bd_edges;   // 1-based indexing (n1, n2, e, l)
    size_t *idx_map;    // 0-based indexing
    double *coords;
    SymBandMatrix *M;
    SymBandMatrix *K;
    SymBandMatrix *M_scalar;
    void (*mesh_model)(double, int);
    void (*set_bk_source)(double, const double[2], double[2]);
    void (*set_bd_disp)(int, char, const double[2], double[1]);
    void (*set_bd_force)(int, char, const double[2], double[2]);
} FE_Model;

typedef struct {
    size_t e1, e2;
    size_t n1, n2;
} Edge;

FE_Model *create_FE_Model(const char *name, ElementType etp, Renumbering rnb);
void free_FE_Model(FE_Model *model);
void load_mesh(FE_Model *model);
void renumber_nodes(FE_Model *model);

void create_edges(
    size_t n_loc,
    size_t n_elem,
    size_t *elems,
    size_t *n_bd_edge_ptr,
    size_t *n_edge_ptr,
    Edge **edges_ptr
);

void find_boundary_nodes(FE_Model *model);

#endif
