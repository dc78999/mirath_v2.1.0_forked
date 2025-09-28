#ifndef MIRATH_TCITH_MIRATH_PARSING_H
#define MIRATH_TCITH_MIRATH_PARSING_H

#include "mirath_parameters.h"
#include "mirath_matrix_ff.h"
#include "mirath_ggm_tree.h"

void unparse_public_key(uint8_t *pk, const seed_t seed_pk, const ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)]);

void parse_public_key(seed_t seed_pk, ff_t y[mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1)], const uint8_t *pk);

void unparse_secret_key(uint8_t *sk, const seed_t seed_sk, const seed_t seed_pk);

void parse_secret_key(seed_t seed_sk, seed_t seed_pk, const uint8_t *sk);


void unparse_signature(uint8_t *signature, const uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t ctr,
                       const hash_t hash2, const mirath_ggm_tree_node_t path[MIRATH_PARAM_MAX_OPEN],
                       const mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                       const ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                       const ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO]);

int parse_signature(uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t *ctr, hash_t hash2,
                    mirath_ggm_tree_node_t path[MIRATH_PARAM_MAX_OPEN], mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                    ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                    ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], const uint8_t *signature);

#endif //MIRATH_TCITH_MIRATH_PARSING_H
