#include <stdint.h>
#include <string.h>

#include "random.h"
#include "mirath_matrix_ff.h"
#include "mirath_parsing.h"
#include "mirath_ggm_tree.h"
#include "mirath_tcith.h"

int mirath_sign(uint8_t *sig_msg, uint8_t *msg, size_t msg_len, uint8_t *sk) {
    uint8_t salt[MIRATH_PARAM_SALT_BYTES] = {0};
    seed_t rseed = {0};
    mirath_ggm_tree_node_t path[MIRATH_PARAM_MAX_OPEN] = {0};
    hash_t h_mpc;

    ff_t S[MIRATH_VAR_FF_S_BYTES] = {0};
    ff_t C[MIRATH_VAR_FF_C_BYTES] = {0};
    ff_t H[MIRATH_VAR_FF_H_BYTES] = {0};

    uint8_t pk[MIRATH_PUBLIC_KEY_BYTES] = {0};

    ff_mu_t S_base[MIRATH_PARAM_TAU][MIRATH_VAR_S] = {0};
    ff_mu_t C_base[MIRATH_PARAM_TAU][MIRATH_VAR_C] = {0};
    ff_mu_t v_base[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO] = {0};
    ff_mu_t v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO] = {0};
    ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES] = {0};
    hash_t h_sh;
    mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU];
    mirath_ggm_tree_t tree = {0};

    ff_mu_t Gamma[MIRATH_VAR_GAMMA];

    ff_mu_t alpha_mid[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO];

    /*
     * Phase 0: Initialization
     * step 1
     */
    mirath_matrix_decompress_secret_key(S, C, H, pk, (const uint8_t *)sk);
    // step 2
    randombytes(salt, MIRATH_PARAM_SALT_BYTES);

    // step 3
    randombytes(rseed, MIRATH_SECURITY_BYTES);

    mirath_ggm_tree_leaves_t seeds = {0};

    memcpy(tree[0], rseed, MIRATH_SECURITY_BYTES);
    mirath_ggm_tree_expand(tree, salt);
    mirath_ggm_tree_get_leaves(seeds, tree); // First output of MultiVC.Commit

    /*
     * Phase 1: Build and Commit to Witness Polynomials
     * step 4
     */
    commit_witness_polynomials(S_base, C_base, v_base, v, h_sh, seeds, aux, salt, S, C);

    /*
     * step 5
     * This block of code refers to Algorithm 9 Challenge matrix Gamma.
     */
    shake_prng_t seedexpander_shake;
    seedexpander_shake_init(&seedexpander_shake, h_sh, 2 * MIRATH_SECURITY_BYTES, NULL, 0);
    seedexpander_shake_get_bytes(&seedexpander_shake, (uint8_t*)Gamma, sizeof(ff_mu_t) * MIRATH_VAR_GAMMA);

    /*
     *Initialization for computing Hash2
     */
    uint8_t domain_separator;
    domain_separator = DOMAIN_SEPARATOR_HASH2_PARTIAL;
    hash_sha3_ctx hash_mpc_ctx;
    hash_init(&hash_mpc_ctx);
    hash_update(&hash_mpc_ctx, &domain_separator, sizeof(uint8_t));
    hash_update(&hash_mpc_ctx, pk, MIRATH_PUBLIC_KEY_BYTES);
    hash_update(&hash_mpc_ctx, salt, MIRATH_PARAM_SALT_BYTES);
    hash_update(&hash_mpc_ctx, msg, msg_len);
    hash_update(&hash_mpc_ctx, h_sh, 2 * MIRATH_SECURITY_BYTES);

    /*
     * Phase 2: MPC simulation
     * steps 6 and 7
     */
    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
        ff_mu_t alpha_base[MIRATH_PARAM_RHO];

        compute_polynomial_proof(alpha_base, alpha_mid[e], S, S_base[e], C, C_base[e], v[e], v_base[e], Gamma, H);

        // Regarding Hash2
        hash_update(&hash_mpc_ctx, (uint8_t*)alpha_base, sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
        hash_update(&hash_mpc_ctx, (uint8_t*)alpha_mid[e], sizeof(ff_mu_t) * MIRATH_PARAM_RHO);
    }

    // step 8
    hash_finalize(h_mpc, &hash_mpc_ctx);

    /*
     * Phase 3: Sharing Opening.
     * step 9
     */
    size_t psi_i_star[MIRATH_PARAM_TAU];

    // SHAKE input (fixed prefix)
    uint8_t shake_input[SHAKE_STEP * (2 * MIRATH_SECURITY_BYTES + sizeof(uint64_t))] = {0};
    for (uint32_t i = 0; i < SHAKE_STEP; i++) {
        memcpy(&shake_input[(2 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)) * i], h_mpc, 2 * MIRATH_SECURITY_BYTES);
    }

    mirath_tcith_view_challenge_t i_star[SHAKE_STEP] = {0};
    uint8_t v_grinding[MIRATH_PARAM_HASH_2_MASK_BYTES * SHAKE_STEP] = {0};

    uint64_t ctr = 0;
    uint32_t valid_step_found = false;

    while(!valid_step_found) {
        // SHAKE input (suffix corresponds with the counter)
        for (uint32_t i = 0; i < SHAKE_STEP; i++) {
            const uint64_t ctr_i = ctr + i;
            memcpy(&shake_input[(2 * MIRATH_SECURITY_BYTES + sizeof(uint64_t)) * i + (2 * MIRATH_SECURITY_BYTES)],
                   (uint8_t *) &ctr_i, sizeof(uint64_t));
        }

        mirath_tcith_expand_view_challenge_x4(i_star, v_grinding, shake_input);

        for (uint32_t step = 0; step < SHAKE_STEP; step++) {
            if (mirath_tcith_discard_input_challenge_2(&v_grinding[step * MIRATH_PARAM_HASH_2_MASK_BYTES])) {
                memset(&v_grinding[step * MIRATH_PARAM_HASH_2_MASK_BYTES], 0, MIRATH_PARAM_HASH_2_MASK_BYTES);
                continue;
            }

            for (size_t e = 0; e < MIRATH_PARAM_TAU; e++) {
                const size_t i = i_star[step][e];
                psi_i_star[e] = mirath_tcith_psi(i, e); // store their respective image under psi
            }

            const uint64_t path_length = mirath_ggm_tree_get_sibling_path(path, (const mirath_ggm_tree_node_t *)tree, psi_i_star);

            if (path_length > MIRATH_PARAM_T_OPEN) {
                memset(path, 0, sizeof(mirath_ggm_tree_node_t) * path_length);
                memset(&v_grinding[step * MIRATH_PARAM_HASH_2_MASK_BYTES], 0, MIRATH_PARAM_HASH_2_MASK_BYTES);
                continue;
            }

            // Valid step found
            ctr += step;
            for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++) {
                const uint32_t idx = mirath_tcith_psi((size_t)i_star[step][e], (size_t)e);
                mirath_commit(&commits_i_star[e], salt, (uint32_t)idx, seeds[idx]);
            }

            valid_step_found = true;
            break;
        }

        if (!valid_step_found) {
            ctr += SHAKE_STEP; // retry with new inputs
        }
    }

    // step 10
    unparse_signature(sig_msg, salt, ctr, h_mpc, path, commits_i_star, aux, alpha_mid);

    return 0;
}
