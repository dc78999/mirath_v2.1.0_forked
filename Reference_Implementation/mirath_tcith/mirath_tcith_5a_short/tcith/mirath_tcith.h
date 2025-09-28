#ifndef MIRATH_MIRATH_TCITH_H
#define MIRATH_MIRATH_TCITH_H

#include <stdint.h>
#include <stdlib.h>

#include "mirath_matrix_ff.h"
#include "mirath_ggm_tree.h"
#include "seed_expand_rijndael.h"
#include "framework_data_types.h"
#include "hash_sha3.h"

#define true 1
#define false 0

static inline size_t mirath_tcith_psi(size_t i, size_t e) {
    if (i < MIRATH_PARAM_N_2) {
        return i * MIRATH_PARAM_TAU + e;
    }
    else {
        return MIRATH_PARAM_N_2 * MIRATH_PARAM_TAU + (i - MIRATH_PARAM_N_2) * MIRATH_PARAM_TAU_1 + e;
    }
}

/* ----------------- Scalar 64-bit chunk ----------------- */
/* little-endian view */
static inline uint64_t load_u64_le(const void *p, size_t avail_bytes) {
    uint64_t v = 0;
    if (avail_bytes >= 8) {
        memcpy(&v, p, 8);
    } else if (avail_bytes > 0) {
        memcpy(&v, p, avail_bytes);
    }
    return v;
}

static inline void store_u64_le(void *p, uint64_t v, size_t avail_bytes) {
    if (avail_bytes >= 8) {
        memcpy(p, &v, 8);
    } else if (avail_bytes > 0) {
        memcpy(p, &v, avail_bytes);
    }
}

/**
* \fn void mirath_tcith_shift_to_left_array(uint8_t *inout_a, size_t length)
* \brief This function performs a shift to right of the input.
*
* \param[in/out] inout_a uint8_t* Representation of a byte string
* \param[in] length size_t Representation of the byte string length
* \param[in] nbits int Representation of the number of bits to be shifted
*/
/* portable 64-bit chunk implementations */
static inline void mirath_tcith_shift_to_right_array(uint8_t *string, const size_t length, const int nbits) {
    if (length == 0 || nbits == 0) return;
    int acc = 0;
    int bits = (nbits > 7) ? 7 : nbits;

    while (acc < nbits) {
        if (length < 16) {
            for (size_t i = 0; i + 1 < length; ++i) {
                string[i] = (uint8_t)((string[i] >> bits) | (string[i + 1] << (8 - bits)));
            }
            string[length - 1] = (uint8_t)(string[length - 1] >> bits);
        } else {
            size_t n_u64 = length / 8;
            size_t tail = length % 8;

            for (size_t i = 0; i + 1 < n_u64; ++i) {
                uint64_t v = load_u64_le(string + i * 8, 8);
                uint64_t next = load_u64_le(string + (i + 1) * 8, 8);
                uint64_t out = (v >> bits) | (next << (64 - bits));
                store_u64_le(string + i * 8, out, 8);
            }
            if (n_u64 > 0) {
                size_t i = n_u64 - 1;
                uint64_t v = load_u64_le(string + i * 8, 8);
                uint64_t next = load_u64_le(string + n_u64 * 8, tail); /* 0..7 bytes */
                uint64_t out = (v >> bits) | (next << (64 - bits));
                store_u64_le(string + i * 8, out, 8);
            }
            size_t start = n_u64 * 8;
            if (start < length) {
                for (size_t j = start; j + 1 < length; ++j) {
                    string[j] = (uint8_t)((string[j] >> bits) | (string[j + 1] << (8 - bits)));
                }
                string[length - 1] = (uint8_t)(string[length - 1] >> bits);
            }
        }
        acc += bits;
        int rem = nbits - acc;
        bits = (rem > 7) ? 7 : rem;
    }
}

static inline void mirath_tcith_shift_to_left_array(uint8_t *string, const int length_in, const int nbits) {
    if (length_in <= 0 || nbits == 0) return;
    size_t length = (size_t)length_in;
    int acc = 0;
    int bits = (nbits > 7) ? 7 : nbits;

    while (acc < nbits) {
        if (length < 16) {
            for (int i = length_in - 1; i >= 1; --i) {
                string[i] = (uint8_t)((string[i] << bits) | (string[i - 1] >> (8 - bits)));
            }
            string[0] = (uint8_t)(string[0] << bits);
        } else {
            size_t n_u64 = length / 8;
            size_t start = n_u64 * 8;

            /* tail bytes (high indices) first */
            if (start < length) {
                for (int j = (int)length - 1; j > (int)start; --j) {
                    string[j] = (uint8_t)((string[j] << bits) | (string[j - 1] >> (8 - bits)));
                }
                if (start > 0) {
                    string[start] = (uint8_t)((string[start] << bits) | (string[start - 1] >> (8 - bits)));
                } else {
                    string[0] = (uint8_t)(string[0] << bits);
                }
            }

            /* full 64-bit words high->low */
            if (n_u64 > 0) {
                for (size_t idx = n_u64; idx-- > 1; ) {
                    size_t i = idx;
                    uint64_t v = load_u64_le(string + i * 8, 8);
                    uint64_t prev = load_u64_le(string + (i - 1) * 8, 8);
                    uint64_t out = (v << bits) | (prev >> (64 - bits));
                    store_u64_le(string + i * 8, out, 8);
                }
                /* first full word */
                {
                    uint64_t v0 = load_u64_le(string + 0, 8);
                    uint64_t out0 = (v0 << bits);
                    store_u64_le(string + 0, out0, 8);
                }
            }
        }
        acc += bits;
        int rem = nbits - acc;
        bits = (rem > 7) ? 7 : rem;
    }
}

/**
* \fn uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *v_grinding)
* \brief This function determines if the w most significant bits of the input are zero.
*
* \param[in] v_grinding String containing the input seed
*/
static inline uint8_t mirath_tcith_discard_input_challenge_2(const uint8_t *v_grinding) {
    uint8_t output = 0x00;
    uint8_t mask = MIRATH_PARAM_HASH_2_MASK;
    for(int i = MIRATH_PARAM_HASH_2_MASK_BYTES - 1; i >= 0 ; i--) {
        output |= (uint8_t)((v_grinding[i] & mask) != 0);
        mask = 0xFF;
    }

    return output;
}

void mirath_commit(mirath_tcith_commit_t *pair_node,
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES],
                   uint32_t i,
                   const uint8_t seed[MIRATH_SECURITY_BYTES]);

void mirath_tcith_internal_steps_pk(ff_t y[MIRATH_VAR_FF_Y_BYTES],
                                    const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                    const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void mirath_tcith_commit(mirath_tcith_commit_t commit, const uint8_t *salt, uint16_t e, uint32_t i, const uint8_t *seed);

void commit_witness_polynomials(ff_mu_t S_base[MIRATH_PARAM_TAU][MIRATH_VAR_S],
                                ff_mu_t C_base[MIRATH_PARAM_TAU][MIRATH_VAR_C],
                                ff_mu_t v_base[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO],
                                ff_mu_t v[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO],
                                hash_t hash_sh,
                                mirath_ggm_tree_leaves_t seeds,
                                ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                                const uint8_t salt[MIRATH_PARAM_SALT_BYTES],
                                const ff_t S[MIRATH_VAR_FF_S_BYTES],
                                const ff_t C[MIRATH_VAR_FF_C_BYTES]);

void compute_share(ff_mu_t S_share[MIRATH_PARAM_TAU][MIRATH_VAR_S], ff_mu_t C_share[MIRATH_PARAM_TAU][MIRATH_VAR_C],
                   ff_mu_t v_share[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], hash_t hash_sh,
                   const mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                   const mirath_tcith_view_challenge_t i_star, const mirath_ggm_tree_leaves_t seeds,
                   const ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                   const uint8_t salt[MIRATH_PARAM_SALT_BYTES]);

void compute_polynomial_proof(ff_mu_t base_alpha[MIRATH_PARAM_RHO], ff_mu_t mid_alpha[MIRATH_PARAM_RHO],
                              const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_mu_t S_rnd[MIRATH_VAR_S],
                              const ff_t C[MIRATH_VAR_FF_C_BYTES], const ff_mu_t C_rnd[MIRATH_VAR_C],
                              const ff_mu_t v[MIRATH_PARAM_RHO], ff_mu_t rnd_v[MIRATH_PARAM_RHO],
                              const ff_mu_t gamma[MIRATH_VAR_GAMMA], const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void recompute_polynomial_proof(ff_mu_t base_alpha[MIRATH_PARAM_RHO], ff_mu_t p,
                                const ff_mu_t S_share[MIRATH_VAR_S], const ff_mu_t C_share[MIRATH_VAR_C],
                                const ff_mu_t v_share[MIRATH_PARAM_RHO], const ff_mu_t gamma[MIRATH_VAR_GAMMA],
                                const ff_t H[MIRATH_VAR_FF_H_BYTES], const ff_t y[MIRATH_VAR_FF_Y_BYTES],
                                const ff_mu_t mid_alpha[MIRATH_PARAM_RHO]);

void mirath_tcith_expand_view_challenge(mirath_tcith_view_challenge_t challenge, uint8_t *v_grinding, const uint8_t *string_input);

#define SHAKE_STEP 4

void mirath_commit_4x(mirath_tcith_commit_t pair_node[SHAKE_STEP],
                      const uint8_t salt[MIRATH_PARAM_SALT_BYTES],
                      const uint32_t i[SHAKE_STEP],
                      const mirath_ggm_tree_node_t seed[SHAKE_STEP]);

#endif //MIRATH_MIRATH_TCITH_H
