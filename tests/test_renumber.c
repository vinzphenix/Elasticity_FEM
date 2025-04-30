/**
 * File:    test_renumber.c
 * Author:  Vincent Degrooff
 * Created: 2025
 * 
 * Description: 
 *   Test the renumbering strategy RCMK
 * 
 * Project: 
 *   FEM Simulation Toolkit for Linear Elasticity
 */

#include "renumber.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// clang-format off
void test_triangle() {
    // TRIANGLE MESH
    size_t n_loc = 3;
    size_t n_elem = 8;
    size_t elems[] = {
        1, 3, 2, 
        2, 3, 4, 
        3, 5, 4, 
        4, 5, 6,
        2, 4, 8, 
        2, 8, 7, 
        4, 6, 9, 
        4, 9, 8
    };
    size_t n_node = 9;
    double coords[] = {
        0.2, 2.,
        1.2, 2.1,
        0.1, 1.,
        1.1, 1.1,
        0., 0.,
        1., 0.1,
        2.2, 2.,
        2.1, 1.,
        2.0, 0.
    };
    size_t n_bd_node = 8;
    size_t bd_nodes[] = {2, 1, 3, 5, 5, 6, 1, 3, 6, 9, 9, 8, 8, 7, 7, 2};
    
    size_t *idx_map;
    renumber(
        n_elem, n_loc, elems, n_bd_node, bd_nodes, n_node, coords, &idx_map, 2
    );
    for (size_t i = 0; i < n_node; i++) {
        printf("%zu -> %zu\n", i + 1, idx_map[i] + 1);
    }
    size_t band = compute_band(n_elem, 3, elems, idx_map);
    printf("Band size = %zu\n", band);
    free(idx_map);
}


void test_quad() {
    // QUAD MESH
    // 1 -- 2 -- 7
    // |  / |  \ |
    // 3 -- 4 -- 8
    // |  / |  \ |
    // 5 -- 6 -- 9
    size_t n_loc = 4;
    size_t n_elem = 3;
    size_t elems[] = {
        1, 2, 3, 4,
        1, 6, 7, 2,
        5, 6, 1, 4
    };
    size_t n_node = 7;
    double coords[] = {
        +0., +0.,
        +1., +1.,
        +0., +2.,
        -1., +1.,
        -1.1, -1.,
        -0.1, -1.,
        +0.9, -1.
    };
    size_t n_bd_node = 6;
    size_t bd_nodes[] = {
        2, 3, 0, 1, 
        3, 4, 0, 2,
        4, 5, 2, 3,
        5, 6, 2, 0,
        6, 7, 1, 1,
        7, 2, 1, 2
    };
    size_t *idx_map;
    renumber(
        n_elem, n_loc, elems, n_bd_node, bd_nodes, n_node, coords, &idx_map, 3
    );
    for (size_t i = 0; i < n_node; i++) {
        printf("%zu -> %zu\n", i + 1, idx_map[i] + 1);
    }
    size_t band = compute_band(n_elem, 3, elems, idx_map);
    printf("Band size = %zu\n", band);
    free(idx_map);
}
// clang-format on

int main(int argc, char *argv[]) {
    test_triangle();
    test_quad();
    return 0;
}
    