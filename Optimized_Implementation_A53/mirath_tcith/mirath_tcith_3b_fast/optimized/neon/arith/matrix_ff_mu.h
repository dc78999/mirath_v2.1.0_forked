#ifndef MIRATH_MATRIX_FF_MU_H
#define MIRATH_MATRIX_FF_MU_H

#include <stdint.h>
#include <string.h>
#include <arm_neon.h>

#include "ff_mu.h"
#include "vector_ff_mu.h"
#include "matrix_ff_arith.h"

#define mirath_matrix_ff_mu_get_entry(m,n,i,j) m[j*n + i]
#define mirath_matrix_ff_mu_set_entry(m,n,i,j,v) m[j*n + i] = v

#define mirath_matrix_ff_mu_bytes_size(x, y) ((x) * (y) * sizeof(ff_mu_t))

#define MIRATH_VAR_FF_MU_S_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_T_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_E_A_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1))
#define MIRATH_VAR_FF_MU_K_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_K, 1))

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

static inline void mirath_matrix_ff_mu_copy(ff_mu_t *matrix1, const ff_mu_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    memcpy(matrix1, matrix2, mirath_matrix_ff_mu_bytes_size(n_rows, n_cols));
}

static inline __attribute__((always_inline)) uint8x16_t read16column(const ff_mu_t *col, const uint32_t bytes) {
    if (bytes == 16) { return vld1q_u8((const uint8_t *)col); }
    uint8_t tmp[16] __attribute__((aligned(16))) = {0};
    for (uint32_t i = 0; i < bytes; ++i) { tmp[i] = col[i]; }
    return vld1q_u8(tmp);
}

static inline __attribute__((always_inline)) void write16column(ff_mu_t *col, const uint8x16_t data, const uint32_t bytes) {
    if (bytes == 16) {
        vst1q_u8((uint8_t *)col, data);
        return;
    }
    uint8_t tmp[16] __attribute__((aligned(16)));
    vst1q_u8(tmp, data);
    for (uint32_t i = 0; i < bytes; ++i) { col[i] = tmp[i]; }
}

static inline void mirath_matrix_ff_mu_add(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                             const ff_mu_t *matrix3, const uint32_t n_rows,
                             const uint32_t n_cols) {
    mirath_vector_ff_mu_add(matrix1, matrix2, matrix3, n_rows*n_cols);
}

static inline void mirath_matrix_ff_mu_add_mu1ff(ff_mu_t *matrix1, const ff_mu_t *matrix2, const ff_t *matrix3,
                                                 const uint32_t n_rows, const uint32_t n_cols) {
    if (n_rows % 8 == 0) {
        mirath_vector_ff_mu_add_ff(matrix1, matrix2, matrix3, n_rows*n_cols);
        return;
    }

    const uint32_t limit = n_rows % 16;

    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = n_rows;
        ff_t *in     = (ff_t    *)matrix3 + j*((n_rows + 1)/2);
        ff_mu_t *in2 = (ff_mu_t *)matrix2 + j*n_rows;
        ff_mu_t *out = matrix1 + j*n_rows;

        while (i >= 16u) {
            uint8x16_t m2 = vld1q_u8((const uint8_t *)in2);
            uint8x16_t t4 = mirath_ff_mu_extend_gf16_x16((const uint8_t *)in);
            uint8x16_t t5 = veorq_u8(t4, m2);
            vst1q_u8((uint8_t *)out, t5);

            in  += 8u;
            out += 16u;
            in2 += 16u;
            i   -= 16u;
        }

        if (limit) {
            uint8_t tmp1[16] __attribute__((aligned(16))) = {0};
            uint8_t tmp2[16] __attribute__((aligned(16))) = {0};
            for (uint32_t k = 0; k < i; k++) { tmp1[k] = in2[k]; }
            uint8x16_t m2 = vld1q_u8(tmp1);

            for (uint32_t k = 0; k < (i + 1) / 2; k++) { tmp2[k] = in[k]; }
            uint8x16_t t4 = mirath_ff_mu_extend_gf16_x16(tmp2);
            uint8x16_t t5 = veorq_u8(t4, m2);
            vst1q_u8(tmp1, t5);
            for (uint32_t k = 0; k < i; k++) { out[k] = tmp1[k]; }
        }
    }
}

static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
                                         const uint32_t n_rows, const uint32_t n_cols) {
    if ((n_rows % 2) == 2) {
        mirath_vector_ff_mu_add_multiple_ff(matrix1, matrix1, scalar, matrix2, n_rows * n_cols);
        return;
    }

    uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(scalar);
    uint8x16_t mask = vdupq_n_u8(0xf);

    const uint32_t limit = n_rows % 16;

    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = n_rows;
        ff_t *in     = (ff_t    *)matrix2 + j*((n_rows + 1)/2);
        ff_mu_t *out = matrix1 + j*n_rows;

        while (i >= 16u) {
            uint8x16_t m2 = vld1q_u8((const uint8_t *)out);
            uint8x16_t t4 = mirath_ff_mu_extend_gf16_x16((const uint8_t *)in);
            uint8x16_t t5 = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], t4, mask);
            uint8x16_t t6 = veorq_u8(t5, m2);
            vst1q_u8((uint8_t *)out, t6);

            in  += 8u;
            out += 16u;
            i   -= 16u;
        }

        if (limit) {
            uint8_t tmp1[16] __attribute__((aligned(16))) = {0};
            uint8_t tmp2[16] __attribute__((aligned(16))) = {0};
            for (uint32_t k = 0; k < limit; k++) { tmp1[k] = out[k]; }
            uint8x16_t m2 = vld1q_u8(tmp1);

            for (uint32_t k = 0; k < (limit + 1) / 2; k++) { tmp2[k] = in[k]; }
            uint8x16_t t4 = mirath_ff_mu_extend_gf16_x16(tmp2);
            uint8x16_t t5 = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], t4, mask);
            uint8x16_t t6 = veorq_u8(t5, m2);
            vst1q_u8(tmp1, t6);
            for (uint32_t k = 0; k < limit; k++) { out[k] = tmp1[k]; }
        }
    }
}

static inline void mirath_matrix_ff_mu_add_multiple_3(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                      const ff_mu_t scalar, const ff_mu_t *matrix3,
                                      const uint32_t n_rows, const uint32_t n_cols) {
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix2, n_rows, i, j);
            const ff_mu_t entry2 = mirath_matrix_ff_mu_get_entry(matrix3, n_rows, i, j);
            const ff_mu_t entry3 = entry1 ^ mirath_ff_mu_mult(scalar, entry2);

            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry3);
        }
    }
}

static inline void mirath_matrix_product_gf16_1_vector_u128(ff_mu_t *result,
                                                        const ff_t *matrix1,
                                                        const ff_mu_t *matrix2,
                                                        const uint32_t n_rows1,
                                                        const uint32_t n_cols1) {
    uint8_t tmp[16] __attribute__((aligned(16)));

    const uint32_t limit = n_rows1 % 16;
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col * bytes_per_col;
        const uint8x16_t b = vdupq_n_u8(*(matrix2 + col));
        ff_mu_t *r = result;

        while ((i + 16) <= n_rows1) {
            uint8x16_t a = mirath_ff_mu_extend_gf16_x16(m1);
            uint8x16_t t = gf256v_mul_u128(b, a);
            t = veorq_u8(t, vld1q_u8((const uint8_t *)r));
            vst1q_u8((uint8_t *)r, t);

            m1 += 8;
            r  += 16;
            i  += 16;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { tmp[j] = m1[j]; }
            uint8x16_t a = mirath_ff_mu_extend_gf16_x16(tmp);
            uint8x16_t t = gf256v_mul_u128(b, a);
            vst1q_u8(tmp, t);
            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp[j]; }
        }
    }
}

static inline void mirath_matrix_ff_mu_product_ff1mu(ff_mu_t *result, const ff_t *matrix1,
                                       const ff_mu_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    if (n_cols2 == 1) {
        memset(result, 0, n_rows1);
        mirath_matrix_product_gf16_1_vector_u128(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_get_entry(matrix1, n_rows1, i, k);
                entry_i_k = mirath_map_ff_to_ff_mu[entry_i_k];
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

static inline void mirath_matrix_ff_mu_product_mu1ff(ff_mu_t *result, const ff_mu_t *matrix1,
                                       const ff_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_get_entry(matrix2, n_cols1, k, j);
                entry_k_j = mirath_map_ff_to_ff_mu[entry_k_j];
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

static inline void mirath_matrix_product_vector_u128(ff_mu_t *result,
                                                  const ff_mu_t *matrix1,
                                                  const ff_mu_t *matrix2,
                                                  const uint32_t n_rows1,
                                                  const uint32_t n_cols1) {
    uint8_t tmp[16] __attribute__((aligned(16)));

    const uint32_t limit = n_rows1 % 16;
    const uint8x16_t m = vdupq_n_u8(0x0F);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const ff_mu_t *m1 = matrix1 + col*n_rows1;

        uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(matrix2[col]);

        ff_mu_t *r = result;
        while ((i + 16) <= n_rows1) {
            uint8x16_t a = vld1q_u8((const uint8_t *)m1);
            uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, m);
            t = veorq_u8(t, vld1q_u8((const uint8_t *)r));
            vst1q_u8((uint8_t *)r, t);

            m1 += 16;
            r  += 16;
            i  += 16;
        }

        if (limit) {
            for (uint32_t j = 0; j < limit; ++j) { tmp[j] = m1[j]; }
            uint8x16_t a = vld1q_u8(tmp);
            uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, m);
            vst1q_u8(tmp, t);
            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp[j]; }
        }
    }
}

static inline void mirath_matrix_product_le32xBxle16_u128(ff_mu_t *C,
                                                         const ff_mu_t *A,
                                                         const ff_mu_t *B,
                                                         const uint32_t a,
                                                         const uint32_t b,
                                                         const uint32_t c) {
    const uint8x16_t m = vdupq_n_u8(0x0F);
    uint8x16_t tmp_C_0[16];
    uint8x16_t tmp_C_1[16];
    for (uint32_t i = 0; i < 16; i++) {
        tmp_C_0[i] = vdupq_n_u8(0);
        tmp_C_1[i] = vdupq_n_u8(0);
    }

    for (uint32_t i = 0; i < b; ++i) {
        uint8x16_t col_A_0 = read16column(A + i*a, MIN(a, 16));
        uint8x16_t col_A_1 = (a > 16) ? read16column(A + i*a + 16, a - 16) : vdupq_n_u8(0);

        for (uint32_t j = 0; j < c; ++j) {
            const ff_mu_t elm_B = B[j*b + i];
            uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(elm_B);

            uint8x16_t r0 = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], col_A_0, m);
            tmp_C_0[j] = veorq_u8(tmp_C_0[j], r0);

            if (a > 16) {
                uint8x16_t r1 = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], col_A_1, m);
                tmp_C_1[j] = veorq_u8(tmp_C_1[j], r1);
            }
        }
    }

    for (uint32_t i = 0; i < c; ++i) {
        write16column(C + i*a, tmp_C_0[i], MIN(a, 16));
        if (a > 16) {
            write16column(C + i*a + 16, tmp_C_1[i], a - 16);
        }
    }
}

static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                 const uint32_t n_rows1, const uint32_t n_cols1,
                                 const uint32_t n_cols2) {
    memset(result, 0, n_rows1*n_cols2);
    if (n_cols2 == 1) {
        mirath_matrix_product_vector_u128(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    } else {
        mirath_matrix_product_le32xBxle16_u128(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
        return;
    }
}

static inline void mirath_matrix_ff_mu_add_product(ff_mu_t *result, const ff_mu_t *matrix1,
                                     const ff_mu_t *matrix2, const uint32_t n_rows1,
                                     const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = mirath_matrix_ff_mu_get_entry(result, n_rows1, i, j);

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

static inline void mirath_matrix_map_ff_to_ff_mu(ff_mu_t *out, const uint8_t *input, const uint32_t nrows, const uint32_t ncols) {
    for (uint32_t i = 0; i < ncols; ++i) {
        for (uint32_t j = 0; j < nrows; ++j) {
            const uint8_t tmp = mirath_matrix_ff_get_entry(input, nrows, j, i);
            *out = mirath_map_ff_to_ff_mu[tmp];
            out += 1;
        }
    }
}

#endif
#undef MIN
