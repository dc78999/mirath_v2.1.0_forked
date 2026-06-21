#include <string.h>
#include "framework_data_types.h"
#include "mirath_arith.h"
#include "mirath_ggm_tree.h"
#include "hash_sha3.h"
#include "mirath_tcith.h"

#define BLOCK_S_SIZE_BYTES ((MIRATH_PARAM_M * MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS * MIRATH_PARAM_TAU + 7) / 8)
#define BLOCK_C_SIZE_BYTES ((MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R) * MIRATH_PARAM_Q_BITS * MIRATH_PARAM_TAU + 7) / 8)
#define BLOCK_ALPHA_SIZE_BYTES ((MIRATH_PARAM_RHO * MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS * MIRATH_PARAM_TAU + 7) / 8)
#define BLOCK_S_SIZE_BITS (MIRATH_PARAM_M * MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS)
#define BLOCK_C_SIZE_BITS (MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R) * MIRATH_PARAM_Q_BITS)
#define BLOCK_ALPHA_SIZE_BITS (MIRATH_PARAM_RHO * MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS)

#define BLOCK_SIZE_BYTES (((BLOCK_S_SIZE_BITS + BLOCK_C_SIZE_BITS + BLOCK_ALPHA_SIZE_BITS) * MIRATH_PARAM_TAU + 7) / 8)

void unparse_public_key(uint8_t *pk, const seed_t seed_pk, const ff_t y[MIRATH_VAR_FF_Y_BYTES]) {
    memcpy(pk, seed_pk, MIRATH_SECURITY_BYTES);

    memcpy(pk + MIRATH_SECURITY_BYTES, y, MIRATH_VAR_FF_Y_BYTES);
}

void parse_public_key(seed_t seed_pk, ff_t y[MIRATH_VAR_FF_Y_BYTES], const uint8_t *pk) {
    memcpy(seed_pk, pk, MIRATH_SECURITY_BYTES);

    memcpy(y, pk + MIRATH_SECURITY_BYTES, MIRATH_VAR_FF_Y_BYTES);
}

void unparse_secret_key(uint8_t *sk, const seed_t seed_sk, const seed_t seed_pk) {
    memcpy(sk, seed_sk, MIRATH_SECURITY_BYTES);

    memcpy(sk + MIRATH_SECURITY_BYTES, seed_pk, MIRATH_SECURITY_BYTES);
}

void parse_secret_key(seed_t seed_sk, seed_t seed_pk, const uint8_t *sk)
{
    memcpy(seed_sk, sk, MIRATH_SECURITY_BYTES);
    memcpy(seed_pk, sk + MIRATH_SECURITY_BYTES, MIRATH_SECURITY_BYTES);
}

void unparse_signature(uint8_t *signature, const uint8_t salt[MIRATH_PARAM_SALT_BYTES], const uint64_t ctr,
                       const hash_sha3_ctx *hash2, const mirath_ggm_tree_node_t path[MIRATH_PARAM_MAX_OPEN],
                       const mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                       const ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                       const ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO])
{
    uint32_t current_byte = 0;

    memcpy(&signature[current_byte], salt, MIRATH_PARAM_SALT_BYTES);
    current_byte += MIRATH_PARAM_SALT_BYTES;

    memcpy(&signature[current_byte], &ctr, sizeof(uint64_t));
    current_byte += sizeof(uint64_t);

    memcpy(&signature[current_byte], hash2, 2 * MIRATH_SECURITY_BYTES);
    current_byte += 2 * MIRATH_SECURITY_BYTES;

    memcpy(&signature[current_byte], path, MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN);
    current_byte += MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN;

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++)
    {
        memcpy(&signature[current_byte], commits_i_star[e], 2 * MIRATH_SECURITY_BYTES);
        current_byte += 2 * MIRATH_SECURITY_BYTES;
    }

    // packing tight field elements
    uint8_t block_s_bytes[BLOCK_S_SIZE_BYTES] = {0};
    uint8_t block_c_bytes[BLOCK_C_SIZE_BYTES] = {0};
    uint8_t block_mu_bytes[BLOCK_ALPHA_SIZE_BYTES] = {0};

    // Pack aux S part
    for (int i = MIRATH_PARAM_TAU - 1; i >= 0; --i)
    {
        const uint32_t block_s_step = mirath_matrix_ff_bytes_per_column(MIRATH_PARAM_M);
        uint32_t block_s_bits = MIRATH_PARAM_M * MIRATH_PARAM_Q_BITS;
        const uint32_t local_s_len = (MIRATH_PARAM_M * MIRATH_PARAM_Q_BITS + 7) / 8;
        const uint32_t local_s_shift = block_s_bits % 8;
        uint8_t local_s_mask = 0xff;
        if (local_s_shift != 0) {
            local_s_mask = (uint8_t)(((uint8_t)1 << local_s_shift) - 1);
        }

        for (int j = MIRATH_PARAM_R - 1; j >= 0; --j)
        {
            mirath_tcith_shift_to_left_array(block_s_bytes, (int)BLOCK_S_SIZE_BYTES, (int)block_s_bits);
            uint8_t local_s_bytes[(MIRATH_PARAM_M * MIRATH_PARAM_Q_BITS + 7) / 8] = {0};
            // We copy bytes of the form: ## .. ## 0#
            memcpy(local_s_bytes, &aux[i][j * block_s_step], local_s_len);
            memcpy(block_s_bytes, local_s_bytes, local_s_len - 1);
            block_s_bytes[local_s_len - 1] ^= (local_s_bytes[local_s_len - 1] & local_s_mask);
        }
    }

    // Pack aux C part
    for (int i = MIRATH_PARAM_TAU - 1; i >= 0; --i)
    {
        const uint32_t block_c_step = mirath_matrix_ff_bytes_per_column(MIRATH_PARAM_R);
        uint32_t block_c_bits = MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS;
        const uint32_t local_c_len = (MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS + 7) / 8;
        const uint32_t local_c_shift = block_c_bits % 8;
        uint8_t local_c_mask = 0xff;
        if (local_c_shift != 0) {
            local_c_mask = (uint8_t)(((uint8_t)1 << local_c_shift) - 1);
        }

        for (int j = (MIRATH_PARAM_N - MIRATH_PARAM_R) - 1; j >= 0; --j)
        {
            mirath_tcith_shift_to_left_array(block_c_bytes, (int)BLOCK_C_SIZE_BYTES, (int)block_c_bits);
            uint8_t local_c_bytes[(MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS + 7) / 8] = {0};
            // We copy bytes of the form: ## .. ## 0#
            memcpy(local_c_bytes, &aux[i][j * block_c_step + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)], local_c_len);
            memcpy(block_c_bytes, local_c_bytes, local_c_len - 1);
            block_c_bytes[local_c_len - 1] ^= (local_c_bytes[local_c_len - 1] & local_c_mask);
        }
    }

    // Pack mid_alpha part
    for (int i = MIRATH_PARAM_TAU - 1; i >= 0; --i)
    {
        const uint32_t local_mu_len = (MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS + 7) / 8;
        uint32_t block_mu_bits = MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS;
        const uint32_t local_mu_shift = block_mu_bits % 8;
        uint8_t local_mu_mask = 0xff;
        if (local_mu_shift != 0) {
            local_mu_mask = (uint8_t)(((uint8_t)1 << local_mu_shift) - 1);
        }
        for (int j = MIRATH_PARAM_RHO - 1; j >= 0; --j)
        {
            mirath_tcith_shift_to_left_array(block_mu_bytes, (int)BLOCK_ALPHA_SIZE_BYTES, (int)block_mu_bits);
            uint8_t local_mu_bytes[(MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS + 7) / 8] = {0};
            // We copy bytes of the form: ## .. ## 0#
            memcpy(local_mu_bytes, (uint8_t *)&mid_alpha[i][j], local_mu_len);
            memcpy(block_mu_bytes, local_mu_bytes, local_mu_len - 1);
            block_mu_bytes[local_mu_len - 1] ^= (local_mu_bytes[local_mu_len - 1] & local_mu_mask);
        }
    }

    memset(&signature[current_byte], 0, BLOCK_SIZE_BYTES);
    memcpy(&signature[current_byte], block_mu_bytes, BLOCK_ALPHA_SIZE_BYTES);
    mirath_tcith_shift_to_left_array(&signature[current_byte], (int)BLOCK_SIZE_BYTES, (int)(BLOCK_C_SIZE_BITS * MIRATH_PARAM_TAU));

    uint8_t mask = 0xff;
    uint32_t shift = (BLOCK_C_SIZE_BITS * MIRATH_PARAM_TAU) % 8;
    if (shift != 0) {
        mask = (uint8_t)(((uint8_t)1 << shift) - 1);
    }

    for (int k = (int)BLOCK_C_SIZE_BYTES - 1; k >= 0; --k) {
        signature[current_byte + k] ^= (block_c_bytes[k] & mask);
        mask = 0xff;
    }

    mirath_tcith_shift_to_left_array(&signature[current_byte], (int)BLOCK_SIZE_BYTES, (int)(BLOCK_S_SIZE_BITS * MIRATH_PARAM_TAU));

    mask = 0xff;
    shift = (BLOCK_S_SIZE_BITS * MIRATH_PARAM_TAU) % 8;
    if (shift != 0) {
        mask = (uint8_t)(((uint8_t)1 << shift) - 1);
    }

    for (int k = (int)BLOCK_S_SIZE_BYTES - 1; k >= 0; --k) {
        signature[current_byte + k] ^= (block_s_bytes[k] & mask);
        mask = 0xff;
    }
}

int parse_signature(uint8_t salt[MIRATH_PARAM_SALT_BYTES], uint64_t *ctr, hash_sha3_ctx *hash2,
                    mirath_ggm_tree_node_t path[MIRATH_PARAM_MAX_OPEN], mirath_tcith_commit_t commits_i_star[MIRATH_PARAM_TAU],
                    ff_t aux[MIRATH_PARAM_TAU][MIRATH_VAR_FF_AUX_BYTES],
                    ff_mu_t mid_alpha[MIRATH_PARAM_TAU][MIRATH_PARAM_RHO], const uint8_t *signature)
{
    // Below code catches trivial forgery
    uint32_t tmp_bits = ((MIRATH_PARAM_M * MIRATH_PARAM_R) + (MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)) + (MIRATH_PARAM_RHO  * MIRATH_PARAM_MU)) * MIRATH_PARAM_TAU;
    tmp_bits *= 4;
    if ((tmp_bits % 8) != 0) {
        const uint8_t mask = (1u << (tmp_bits % 8)) - 1;
        if ((signature[MIRATH_SIGNATURE_BYTES - 1] & mask) != signature[MIRATH_SIGNATURE_BYTES - 1]) {
            return 1;
        }
    }

    uint32_t current_byte = 0;

    memcpy(salt, &signature[current_byte], MIRATH_PARAM_SALT_BYTES);
    current_byte += MIRATH_PARAM_SALT_BYTES;

    memcpy(ctr, &signature[current_byte], sizeof(uint64_t));
    current_byte += sizeof(uint64_t);

    memcpy(hash2, &signature[current_byte], 2 * MIRATH_SECURITY_BYTES);
    current_byte += 2 * MIRATH_SECURITY_BYTES;

    memcpy(path, &signature[current_byte], MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN);
    current_byte += MIRATH_SECURITY_BYTES * MIRATH_PARAM_T_OPEN;

    for (uint32_t e = 0; e < MIRATH_PARAM_TAU; e++)
    {
        memcpy(commits_i_star[e], &signature[current_byte], 2 * MIRATH_SECURITY_BYTES);
        current_byte += 2 * MIRATH_SECURITY_BYTES;
    }

    // Calculate block parameters
    uint8_t block_bytes[BLOCK_SIZE_BYTES] = {0};
    uint8_t block_s_bytes[BLOCK_S_SIZE_BYTES] = {0};
    uint8_t block_c_bytes[BLOCK_C_SIZE_BYTES] = {0};

    uint8_t block_mu_bytes[BLOCK_ALPHA_SIZE_BYTES] = {0};

    memcpy(block_bytes, &signature[current_byte], BLOCK_SIZE_BYTES);
    memcpy(block_s_bytes, block_bytes, BLOCK_S_SIZE_BYTES);

    // Pack aux S part
    for (size_t i = 0; i < MIRATH_PARAM_TAU; i++)
    {
        const uint32_t block_s_step = mirath_matrix_ff_bytes_per_column(MIRATH_PARAM_M);
        const uint32_t local_s_len = (MIRATH_PARAM_M * MIRATH_PARAM_Q_BITS + 7) / 8;
        uint32_t block_s_bits = MIRATH_PARAM_M * MIRATH_PARAM_Q_BITS;
        const uint32_t local_s_shift = block_s_bits % 8;
        uint8_t local_s_mask = 0xff;
        if (local_s_shift != 0) {
            local_s_mask = (uint8_t)(((uint8_t)1 << local_s_shift) - 1);
        }

        for (int j = 0; j < MIRATH_PARAM_R; ++j)
        {
            // We copy bytes of the form: ## .. ## 0#
            memcpy(&aux[i][j * block_s_step], block_s_bytes, local_s_len);
            aux[i][j * block_s_step + local_s_len - 1] &= local_s_mask;
            mirath_tcith_shift_to_right_array(block_s_bytes, (int)BLOCK_S_SIZE_BYTES, (int)block_s_bits);
        }
    }

    mirath_tcith_shift_to_right_array(block_bytes, (int)BLOCK_SIZE_BYTES, (int)(BLOCK_S_SIZE_BITS * MIRATH_PARAM_TAU));;
    memcpy(block_c_bytes, block_bytes, BLOCK_C_SIZE_BYTES);

    // Pack aux C part
    for (int i = 0; i < MIRATH_PARAM_TAU; ++i)
    {
        const uint32_t block_c_step = mirath_matrix_ff_bytes_per_column(MIRATH_PARAM_R);
        const uint32_t local_c_len = (MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS + 7) / 8;
        uint32_t block_c_bits = MIRATH_PARAM_R * MIRATH_PARAM_Q_BITS;
        const uint32_t local_c_shift = block_c_bits % 8;
        uint8_t local_c_mask = 0xff;
        if (local_c_shift != 0) {
            local_c_mask = (uint8_t)(((uint8_t)1 << local_c_shift) - 1);
        }

        for (int j = 0; j < (MIRATH_PARAM_N - MIRATH_PARAM_R); ++j)
        {
            // We copy bytes of the form: ## .. ## 0#
            memcpy(&aux[i][j * block_c_step + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R)], block_c_bytes, local_c_len);
            aux[i][j * block_c_step + mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + local_c_len- 1] &= local_c_mask;
            mirath_tcith_shift_to_right_array(block_c_bytes, (int)BLOCK_C_SIZE_BYTES, (int)block_c_bits);
        }
    }

    mirath_tcith_shift_to_right_array(block_bytes, (int)BLOCK_SIZE_BYTES, (int)(BLOCK_C_SIZE_BITS * MIRATH_PARAM_TAU));;
    memcpy(block_mu_bytes, block_bytes, BLOCK_ALPHA_SIZE_BYTES);

    // Pack mid_alpha part
    for (int i = 0; i < MIRATH_PARAM_TAU; ++i)
    {
        uint32_t block_mu_bits = MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS;
        const uint32_t local_mu_len = (MIRATH_PARAM_MU * MIRATH_PARAM_Q_BITS + 7) / 8;
        const uint32_t local_mu_shift = block_mu_bits % 8;
        uint8_t local_mu_mask = 0xff;
        if (local_mu_shift != 0) {
            local_mu_mask = (uint8_t)(((uint8_t)1 << local_mu_shift) - 1);
        }

        for (int j = 0; j < MIRATH_PARAM_RHO; ++j)
        {
            // We copy bytes of the form: ## .. ## 0#
            memcpy((uint8_t *)&mid_alpha[i][j], block_mu_bytes, local_mu_len);
            ((uint8_t *)&mid_alpha[i][j])[local_mu_len - 1] &= local_mu_mask;
            mirath_tcith_shift_to_right_array(block_mu_bytes, (int)BLOCK_ALPHA_SIZE_BYTES, (int)block_mu_bits);
        }
    }
    return 0;
}
