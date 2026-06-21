#ifndef MATRIX_FF_ARITH_H
#define MATRIX_FF_ARITH_H

#include <stdint.h>
#include <arm_neon.h>
#include "ff.h"

#define mirath_matrix_ff_bytes_per_column(n_rows) (((n_rows) >> 1u) + ((n_rows) & 1u))
#define mirath_matrix_ff_bytes_size(n_rows, n_cols) ((mirath_matrix_ff_bytes_per_column(n_rows)) * (n_cols))

#define MIRATH_VAR_FF_AUX_BITS ((MIRATH_PARAM_M * MIRATH_PARAM_R + MIRATH_PARAM_R * (MIRATH_PARAM_N - MIRATH_PARAM_R)) * MIRATH_PARAM_Q_BITS)
#define MIRATH_VAR_FF_AUX_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_S_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R))
#define MIRATH_VAR_FF_C_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_H_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K))
#define MIRATH_VAR_FF_Y_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1))
#define MIRATH_VAR_FF_T_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_E_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N, 1))

#define OFF_E_A ((8 * MIRATH_VAR_FF_Y_BYTES) - 4 * (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K))
#define OFF_E_B ((8 * mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1)) - 4 * MIRATH_PARAM_K)

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

static inline uint32_t gf16_matrix_load4(const ff_t *m,
                           const uint32_t nrows,
                           const uint32_t i,
                           const uint32_t j,
                           const uint32_t bytes) {
    uint8_t tmp[4] = {0};
    const uint32_t pos = mirath_matrix_ff_bytes_per_column(nrows) * j + i/2;
    for (uint32_t k = 0; k < MIN(bytes, 4); k++) {
        tmp[k] = m[pos + k];
    }
    return *((uint32_t*)tmp);
}

static inline void mirath_matrix_set_to_ff(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    if (n_rows & 1) {
        const uint32_t matrix_height =  mirath_matrix_ff_bytes_per_column(n_rows);
        const uint32_t matrix_height_x = matrix_height -  1;
        for (uint32_t i = 0; i < n_cols; i++) {
            matrix[i * matrix_height + matrix_height_x ] &= 0x0f;
        }
    }
}

static inline ff_t mirath_matrix_ff_get_entry(const ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    if (i & 1u) { return  matrix[nbytes_col * j + (i >> 1)] >> 4; }
    else        { return matrix[nbytes_col * j + (i >> 1)] & 0x0f; }
}

static inline void mirath_matrix_ff_set_entry(ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j, const ff_t scalar) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    const uint32_t target_byte_id = nbytes_col * j + (i >> 1);
    if (i & 1) {
        matrix[target_byte_id] &= 0x0f;
        matrix[target_byte_id] |= (scalar << 4);
    } else {
        matrix[target_byte_id] &= 0xf0;
        matrix[target_byte_id] |= scalar;
    }
}

static inline void mirath_matrix_ff_add_arith(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
        const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    for (uint32_t i = 0; i < n_bytes; i++) {
        matrix1[i] = matrix2[i] ^ matrix3[i];
    }
}

static inline void mirath_matrix_ff_add_multiple_arith(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    for (uint32_t i = 0; i < n_bytes; i++) {
        matrix1[i] ^= mirath_ff_product(scalar, matrix2[i] & 0xf);
        matrix1[i] ^= mirath_ff_product(scalar, matrix2[i] >> 4) << 4;
    }
}

static inline void mirath_matrix_ff_product_vector_u128(ff_t *result,
                                                        const ff_t *matrix1,
                                                        const ff_t *matrix2,
                                                        const uint32_t n_rows1,
                                                        const uint32_t n_cols1) {
    uint8_t tmp[16] __attribute__((aligned(16))) = {0};
    const uint32_t limit = n_rows1 % 32;
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);
    const uint8x16_t mask = vdupq_n_u8(0xf);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col * bytes_per_col;
        const uint8_t b2 = mirath_matrix_ff_get_entry(matrix2, 0, col, 0);
        const uint8x16_t bl = mirath_ff_tbl16_multab(b2);
        const uint8x16_t bh = vshlq_n_u8(bl, 4);
        ff_t *r = result;

        while ((i + 32) <= n_rows1) {
            uint8x16_t a = vld1q_u8((const uint8_t *)m1);
            uint8x16_t t = mirath_ff_linear_transform_8x8_128b(bl, bh, a, mask);

            t = veorq_u8(t, vld1q_u8((const uint8_t *)r));
            vst1q_u8((uint8_t *)r, t);

            m1 += 16;
            r  += 16;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { tmp[j] = m1[j]; }
            uint8x16_t a = vld1q_u8(tmp);
            uint8x16_t t = mirath_ff_linear_transform_8x8_128b(bl, bh, a, mask);
            vst1q_u8(tmp, t);
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { r[j] ^= tmp[j]; }
        }
    }
}

static inline void mirath_matrix_ff_product_arith(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {

    if (n_cols2 == 1) {
        mirath_matrix_ff_product_vector_u128(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }

    uint32_t i, j, k;
    ff_t entry_i_k, entry_k_j, entry_i_j;

    const uint32_t matrix_height =  mirath_matrix_ff_bytes_per_column(n_rows1);
    const uint32_t matrix_height_x = matrix_height -  1;

    for (i = 0; i < n_rows1; i++) {
        for (j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_product(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }

    if (n_rows1 & 1) {
        for (i = 0; i < n_cols2; i++) {
            result[i * matrix_height + matrix_height_x] &= 0x0f;
        }
    }
}

#undef MIN
#endif
