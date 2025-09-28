/**
 * @file mirath_ggm_tree.h
 * @brief Header file for mirath_ggm_tree.c
 */

#ifndef MIRATH_GGM_TREE_H
#define MIRATH_GGM_TREE_H

#include <stdio.h>
#include <stdint.h>
#include "mirath_parameters.h"

#define MIRATH_LEAVES_SEEDS_OFFSET (MIRATH_PARAM_TREE_LEAVES - 1)

// No control on the path length for random instances, but experiments suggest 3 * MIRATH_PARAM_T_OPEN as upper bound
#define MIRATH_PARAM_MAX_OPEN (3 * MIRATH_PARAM_T_OPEN)

#if (MIRATH_PARAM_TREE_LEAVES > 0xFFFFFFFF)
#error MIRATH_PARAM_TREE_LEAVES must fit in uint32_t
#endif

typedef uint8_t mirath_ggm_tree_node_t[MIRATH_SECURITY_BYTES];
typedef mirath_ggm_tree_node_t mirath_ggm_tree_t[2 * MIRATH_PARAM_TREE_LEAVES - 1] __attribute__((aligned(16)));
typedef mirath_ggm_tree_node_t mirath_ggm_tree_leaves_t[MIRATH_PARAM_TREE_LEAVES];

void mirath_ggm_tree_expand(mirath_ggm_tree_t ggm_tree, const uint8_t salt[MIRATH_PARAM_SALT_BYTES]);

int mirath_ggm_tree_partial_expand(mirath_ggm_tree_t partial_ggm_tree,
                                 const uint8_t salt[MIRATH_PARAM_SALT_BYTES],
                                 const mirath_ggm_tree_node_t path_seeds[MIRATH_PARAM_MAX_OPEN],
                                 size_t path_length,
                                 const size_t hidden_leaves[MIRATH_PARAM_TAU]);

int mirath_ggm_tree_get_sibling_path(mirath_ggm_tree_node_t path_seeds[MIRATH_PARAM_MAX_OPEN],
                                   const mirath_ggm_tree_t ggm_tree,
                                   const size_t hidden_leaves[MIRATH_PARAM_TAU]);

void mirath_ggm_tree_get_leaves(mirath_ggm_tree_leaves_t output, mirath_ggm_tree_t tree);

void mirath_ggm_tree_print_sibling_path(const mirath_ggm_tree_node_t path[MIRATH_PARAM_T_OPEN]);

#endif //MIRATH_GGM_TREE_H
