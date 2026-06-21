#ifndef MIRATH_MATRIX_FF_MU_12_H
#define MIRATH_MATRIX_FF_MU_12_H

#include <stdint.h>

#include "ff_mu.h"
#include "matrix_ff_arith.h"
#include "vector_ff_mu.h"


#define mirath_matrix_ff_mu_get_entry(m,n,i,j) m[j*n + i]
#define mirath_matrix_ff_mu_set_entry(m,n,i,j,v) m[j*n + i] = v

#define mirath_matrix_ff_mu_bytes_size(x, y) ((x) * (y) * sizeof(ff_mu_t))

#define MIRATH_VAR_FF_MU_S_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_T_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_MU_E_A_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1))
#define MIRATH_VAR_FF_MU_K_BYTES (mirath_matrix_ff_mu_bytes_size(MIRATH_PARAM_K, 1))

/// matrix1 = matrix2
static inline void mirath_matrix_ff_mu_copy(ff_mu_t *matrix1, const ff_mu_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    memcpy(matrix1, matrix2, mirath_matrix_ff_mu_bytes_size(n_rows, n_cols));
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add(ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                           const ff_mu_t *matrix3, const uint32_t n_rows,
 *                           const uint32_t n_cols)
 * \brief matrix1 = matrix2 + matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] matrix3 Matrix over ff_mu
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                             const ff_mu_t *matrix3, const uint32_t n_rows,
                             const uint32_t n_cols) {
    mirath_vector_ff_mu_add(matrix1, matrix2, matrix3, n_rows*n_cols);
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_mu1ff(ff_mu_t *matrix1, const ff_mu_t *matrix2, const ff_t *matrix3,
 *                                               const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 = matrix2 + matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] matrix3 Matrix over ff
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_mu1ff(ff_mu_t *matrix1, const ff_mu_t *matrix2, const ff_t *matrix3,
                                                 const uint32_t n_rows, const uint32_t n_cols) {
    if (n_rows % 8 == 0) {
        mirath_vector_ff_mu_add_ff(matrix1, matrix2, matrix3, n_cols*n_rows);
        return;
    }
    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(n_rows);
    uint16_t tmp[8] __attribute__((aligned(16))) = {0};
    for (uint32_t col = 0; col < n_cols; ++col) {
        uint32_t i = n_rows;
        const ff_t *in2 = matrix3 + col*gf2_col_bytes;
        const ff_mu_t *in1 = matrix2 + col*n_rows;
        ff_mu_t *out = matrix1 + col*n_rows;
        while (i >= 8u) {
            const __m128i m2 = mirath_ff_mu_expand_ff_x8_u256(in2);
            const __m128i m1 = _mm_loadu_si128((__m128i *)in1);
            _mm_storeu_si128((__m128i *)out, m1 ^ m2);
            i   -= 8u;
            in2 += 1u;
            in1 += 8u;
            out += 8u;
        }
        if (i) {
            for (uint32_t j = 0; j < i; j++) { tmp[j] = in1[j];}
            const __m128i m1 = _mm_loadu_si128((__m128i *)tmp);
            const __m128i m2 = mirath_ff_mu_expand_ff_x8_u256(in2);
            _mm_storeu_si128((__m128i *)tmp, m2 ^ m1);
            for (uint32_t j = 0; j < i; j++) { out[j] = tmp[j];}
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
 *                                       const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 += scalar * matrix2
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] scalar scalar over ff_mu
 * \param[in] matrix2 Matrix over ff
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
                                         const uint32_t n_rows, const uint32_t n_cols) {
    if ((n_rows % 8 == 0) || n_cols == 1) {
        mirath_vector_ff_mu_add_multiple_ff(matrix1, matrix1, scalar, matrix2, n_rows*n_cols);
        return;
    }

    const __m256i s256 = _mm256_set1_epi16(scalar);
    const __m128i s128 = _mm_set1_epi16(scalar);
    if (n_rows == 4) {
        const uint64_t m = 0x01010101;
        const uint32_t limit = n_cols % 4;

        const ff_t *mm2 = matrix2;
        ff_mu_t *out = matrix1;

        /// NOTE: asumes n_cols %2 == 0
        uint32_t col = 0;
        for (; (col+4) <= n_cols; col+=4) {
            const uint32_t a = *(uint32_t *)(matrix2 + col);

            const uint64_t t41 = _pdep_u64(a>> 0u, m);
            const uint64_t t42 = _pdep_u64(a>> 8u, m);
            const uint64_t t43 = _pdep_u64(a>>16u, m);
            const uint64_t t44 = _pdep_u64(a>>24u, m);
            const uint64_t t31 =  t41 ^ (t42 << 32);
            const uint64_t t32 =  t43 ^ (t44 << 32);
            const __m128i t21  = _mm_set_epi64x(t32, t31);
            const __m256i t1   = _mm256_cvtepu8_epi16(t21);

            const __m256i m1 = _mm256_loadu_si256((__m256i *)out);
            const __m256i c1 = mirath_ff_mu_mul_ff_u256(s256, t1);
            _mm256_storeu_si256((__m256i *)(out +  0), m1 ^ c1);

            out += 4*n_rows;
            mm2 += 4;
        }

        if (limit) {
            ff_mu_t tmp[16] __attribute__((aligned(32)));

            const uint64_t t41 = _pdep_u64(mm2[0], m);
            const uint64_t t42 = _pdep_u64(mm2[1], m);
            const uint64_t t31 = t41 ^ (t42 << 32);
            const __m128i t21  = _mm_set_epi64x(0, t31);
            const __m128i t1   = _mm_cvtepu8_epi16(t21);
            const __m128i c1 = mirath_ff_mu_mul_ff_u128(s128, t1);
            _mm_store_si128((__m128i *)tmp, c1);

            for (uint32_t i = 0; i < 2*n_rows; i++) {
                out[i] ^= tmp[i];
            }
        }
        return;
    }

    if (n_rows == 5) {
        const uint64_t m = 0x0101010101ull;
        const uint64_t mk = m * ((1ull<<n_rows) - 1ull);
        const uint32_t limit = n_cols % 4u;
        const uint64_t ss = scalar;

        const ff_t *mm2 = matrix2;
        ff_mu_t *out = matrix1;

        /// NOTE: asumes n_cols % 4 == 1
        uint32_t col = 0;
        for (; (col+4) <= n_cols; col+=4) {
            const uint64_t a = *(uint32_t *)(mm2);
            const uint64_t b = _pext_u64(a, mk);

            const uint64_t t21 = _pdep_u64(b     , 0x0101010101010101ull);
            const uint64_t t22 = _pdep_u64(b>> 8u, 0x0101010101010101ull);
            const uint64_t t23 = _pdep_u64(b>>16u, 0x0001000100010001ull);
            const uint64_t s64 = t23 * ss;

            const __m128i t2 = _mm_set_epi64x(t22, t21);
            const __m256i t1 = _mm256_cvtepu8_epi16(t2);
            const __m256i m1 = _mm256_loadu_si256((__m256i *)out);

            const __m256i c1 = mirath_ff_mu_mul_ff_u256(s256, t1);

            _mm256_storeu_si256((__m256i *)out, m1 ^ c1);
            *(uint64_t *)(out + 16) ^= s64;

            out += 4*n_rows;
            mm2 += 4;
        }

        /// NOTE: assumes limit == 1
        if (limit) {
            ff_mu_t tmp[16] __attribute__((aligned(32)));

            const uint64_t a = *(uint8_t *)(mm2);
            const uint64_t b = _pext_u64(a, mk);

            const uint64_t t21 = _pdep_u64(b     , 0x0101010101010101ull);
            // const uint64_t t22 = _pdep_u64(b>> 8u, 0x0101010101010101);

            const __m128i t2 = _mm_set_epi64x(0, t21);
            const __m256i t1 = _mm256_cvtepu8_epi16(t2);
            const __m256i c1 = mirath_ff_mu_mul_ff_u256(s256, t1);

            _mm256_store_si256((__m256i *)tmp, c1);

            for (uint32_t i = 0; i < limit*n_rows; i++) {
                out[i] ^= tmp[i];
            }
        }

        return;
    }

    if (n_rows == 6) {
        const uint64_t m = 0x010101010101;
        const uint32_t limit = n_cols % 4;

        const ff_t *mm2 = matrix2;
        ff_mu_t *out = matrix1;

        /// NOTE: asumes n_cols %2 == 0
        uint32_t col = 0;
        for (; (col+4) <= n_cols; col+=4) {
            const uint32_t a = *(uint32_t *)(mm2);

            const uint64_t t41 = _pdep_u64(a>> 0u, m);
            const uint64_t t42 = _pdep_u64(a>> 8u, m);
            const uint64_t t43 = _pdep_u64(a>>16u, m);
            const uint64_t t44 = _pdep_u64(a>>24u, m);
            const uint64_t t31 = t41 ^ (t42 << 48);
            const uint64_t t32 = (t42 >> 16) ^ (t43 << 32);
            const uint64_t t33 = (t43 >> 32) ^ (t44 << 16);
            const __m128i t21  = _mm_set_epi64x(t32, t31);
            const __m128i t22  = _mm_set_epi64x(0, t33);
            const __m256i t1   = _mm256_cvtepu8_epi16(t21);
            const __m128i t2   = _mm_cvtepi8_epi16(t22);

            const __m256i m1 = _mm256_loadu_si256((__m256i *)out);
            const __m128i m2 = _mm_loadu_si128((__m128i *)(out + 16));

            const __m256i c1 = mirath_ff_mu_mul_ff_u256(s256, t1);
            const __m128i c2 = mirath_ff_mu_mul_ff_u128(s128, t2);

            _mm256_storeu_si256((__m256i *)(out +  0), m1 ^ c1);
            _mm_storeu_si128((__m128i *)(out + 16), m2 ^ c2);

            out += 4*n_rows;
            mm2 += 4;
        }

        if (limit) {
            ff_mu_t tmp[16] __attribute__((aligned(32)));

            const uint64_t t41 = _pdep_u64(mm2[0], m);
            const uint64_t t42 = _pdep_u64(mm2[1], m);
            const uint64_t t31 = t41 ^ (t42 << 48);
            const uint64_t t32 = t42 >> 16;
            const __m128i t21  = _mm_set_epi64x(t32, t31);
            const __m256i t1   = _mm256_cvtepu8_epi16(t21);
            const __m256i c1 = mirath_ff_mu_mul_ff_u256(s256, t1);
            _mm256_store_si256((__m256i *)tmp, c1);

            for (uint32_t i = 0; i < 2*n_rows; i++) {
                out[i] ^= tmp[i];
            }
        }
        return;
    }


    const uint32_t limit = n_rows % 8;
    const __m128i s = _mm_set1_epi16(scalar);
    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(n_rows);
    uint16_t tmp[8] __attribute__((aligned(32))) = {0};

    for (uint32_t col = 0; col < n_cols; ++col) {
        uint32_t i = n_rows;
        const ff_t *in2 = matrix2 + col*gf2_col_bytes;
        ff_mu_t *out = matrix1 + col*n_rows;

        while (i >= 16u) {
            const __m256i m1 = _mm256_loadu_si256((__m256i *)out);
            const __m256i m2 = mirath_ff_mu_expand_ff_x16_u256(in2);
            const __m256i t = mirath_ff_mu_mul_ff_u256(s256, m2);
            _mm256_storeu_si256((__m256i *)out, m1 ^ t);
            i   -= 16u;
            in2 += 2u;
            out += 16u;
        }

        while (i >= 8u) {
            const __m128i m1 = _mm_loadu_si128((__m128i *)out);
            const __m128i m2 = mirath_ff_mu_expand_ff_x8_u256(in2);
            const __m128i t = mirath_ff_mu_mul_ff_u128(s, m2);
            _mm_storeu_si128((__m128i *)out, m1 ^ t);
            i   -= 8u;
            in2 += 1u;
            out += 8u;
        }

        if (limit) {
            for (uint32_t j = 0; j < i; j++) { tmp[j] = out[j];}
            const __m128i m1 = _mm_loadu_si128((__m128i *)tmp);
            const __m128i m2 = mirath_ff_mu_expand_ff_x8_u256(in2);
            const __m128i t = mirath_ff_mu_mul_ff_u128(s, m2);
            _mm_storeu_si128((__m128i *)tmp, m1 ^ t);
            for (uint32_t j = 0; j < i; j++) { out[j] = tmp[j];}
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_multiple_3(ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                                    const ff_mu_t scalar, const ff_mu_t *matrix3,
 *                                    const uint32_t n_rows, const uint32_t n_cols)
 * \brief matrix1 = matrix2 + scalar * matrix3
 *
 * \param[out] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] scalar scalar over ff_mu
 * \param[in] matrix3 Matrix over ff_mu
 * \param[in] n_rows number of rows
 * \param[in] n_cols number of columns
 */
static inline void mirath_matrix_ff_mu_add_multiple_3(ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                      const ff_mu_t scalar, const ff_mu_t *matrix3,
                                      const uint32_t n_rows, const uint32_t n_cols) {
    mirath_vector_ff_mu_add_multiple(matrix1, matrix2, scalar, matrix3, n_rows*n_cols);
}

/// \brief result = matrix1 * matrix2
/// NOTE: assumes that n_cols2 == 1
/// \param[out] result Matrix over ff_mu
/// \param[in] matrix1 Matrix over ff
/// \param[in] matrix2 Matrix over ff_mu
/// \param[in] n_rows1 number of rows in matrix1
/// \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
static inline void mirath_matrix_mul_ff_mu_vector_u256(ff_mu_t *result,
                                                       const ff_t *matrix1,
                                                       const ff_mu_t *matrix2,
                                                       const uint32_t n_rows1,
                                                       const uint32_t n_cols1) {
    const __m256i mask = _mm256_setzero_si256();
    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(n_rows1);
    const uint32_t limit = n_rows1 % 32;

    uint8_t tmp[64] __attribute__((aligned(32)));
    uint16_t *tmp1 = (uint16_t *)tmp;

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col*gf2_col_bytes;
        const __m256i b = _mm256_set1_epi16(*(matrix2 + col));
        ff_mu_t *r = result;
        while ((i + 32) <= n_rows1) {
            const __m256i a1 = mirath_ff_mu_expand_ff_x16_u256(m1 + 0);
            const __m256i a2 = mirath_ff_mu_expand_ff_x16_u256(m1 + 2);

            const __m256i t1 = _mm256_cmpgt_epi16(a1, mask);
            const __m256i t2 = _mm256_cmpgt_epi16(a2, mask);

            const __m256i c1 = _mm256_loadu_si256((const __m256i *)(r +  0));
            const __m256i c2 = _mm256_loadu_si256((const __m256i *)(r + 16));

            const __m256i d1 = c1 ^ (b & t1);
            const __m256i d2 = c2 ^ (b & t2);

            _mm256_storeu_si256((__m256i *)(r +  0), d1);
            _mm256_storeu_si256((__m256i *)(r + 16), d2);
            m1 += 4;
            r  += 32;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+7)/8; j++) { tmp[j] = m1[j]; }
            const __m256i a1 = mirath_ff_mu_expand_ff_x16_u256(tmp + 0);
            const __m256i a2 = mirath_ff_mu_expand_ff_x16_u256(tmp + 2);

            const __m256i t1 = _mm256_cmpgt_epi16(a1, mask);
            const __m256i t2 = _mm256_cmpgt_epi16(a2, mask);

            const __m256i d1 = b & t1;
            const __m256i d2 = b & t2;

            _mm256_storeu_si256((__m256i *)(tmp1 +  0), d1);
            _mm256_storeu_si256((__m256i *)(tmp1 + 16), d2);
            for (uint32_t j = 0; j < limit; j++) { r[j] ^= tmp1[j]; }

        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product_ff1mu(ff_mu_t *result, const ff_t *matrix1,
 *                                     const ff_mu_t *matrix2, const uint32_t n_rows1,
 *                                     const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product_ff1mu(ff_mu_t *result, const ff_t *matrix1,
                                       const ff_mu_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    if (n_cols2 == 1) {
        memset(result, 0, mirath_matrix_ff_mu_bytes_size(n_rows1, n_cols2));
        mirath_matrix_mul_ff_mu_vector_u256(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_get_entry(matrix1, n_rows1, i, k) & 0x01;
                entry_k_j = mirath_matrix_ff_mu_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product_mu1ff(ff_mu_t *result, const ff_mu_t *matrix1,
 *                                     const ff_t *matrix2, const uint32_t n_rows1,
 *                                     const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product_mu1ff(ff_mu_t *result, const ff_mu_t *matrix1,
                                       const ff_t *matrix2, const uint32_t n_rows1,
                                       const uint32_t n_cols1, const uint32_t n_cols2) {
    ff_mu_t entry_i_k, entry_k_j, entry_i_j;

    for(uint32_t i = 0; i < n_rows1; i++) {
        for (uint32_t j = 0; j < n_cols2; j++) {
            entry_i_j = 0;

            for (uint32_t k = 0; k < n_cols1; k++) {
                entry_i_k = mirath_matrix_ff_mu_get_entry(matrix1, n_rows1, i, k);
                entry_k_j = mirath_matrix_ff_get_entry(matrix2, n_cols1, k, j) & 0x01;
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }
}

/// \brief result = matrix1 * matrix2
/// NOTE: assumes `n_cols2` = 1
/// \param[out] result Matrix over ff_mu
/// \param[in] matrix1 Matrix over ff_mu
/// \param[in] matrix2 Matrix over ff_mu
/// \param[in] n_rows1 number of rows in matrix1
/// \param[in] n_cols1 number of columns and rows in matrix1 and matrix2
///                     respectively
static inline void mirath_matrix_ff_mu_mul_u256_vector(ff_mu_t *result,
                                                       const ff_mu_t *matrix1,
                                                       const ff_mu_t *matrix2,
                                                       const uint32_t n_rows1,
                                                       const uint32_t n_cols1) {
    uint16_t tmp[16] __attribute__((aligned(32)));

    const uint32_t limit = n_rows1 % 16;
    const __m256i mask = _mm256_set1_epi16(0xFFF);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint16_t *m1 = (uint16_t *)(matrix1 + col*n_rows1);
        const __m256i b = _mm256_set1_epi16(*(matrix2 + col));
        ff_mu_t *r = result;
        while ((i + 16) <= n_rows1) {
            __m256i a = _mm256_loadu_si256((const __m256i *)m1);

            a &= mask;
            __m256i t = mirath_ff_mu_mul_u256(a, b);
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 16;
            r  += 16;
            i  += 16;
        }

        // tail mngt
        if (limit) {
            for (uint32_t j = 0; j < limit; ++j) { tmp[j] = m1[j]; }
            __m256i a = _mm256_loadu_si256((const __m256i *)tmp);
            a &= mask;
            __m256i t = mirath_ff_mu_mul_u256(a, b);

            _mm256_storeu_si256((__m256i *)tmp, t);

            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp[j]; }
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
 *                               const uint32_t n_rows1, const uint32_t n_cols1,
 *                               const uint32_t n_cols2)
 * \brief result = matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_add_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                 const uint32_t n_rows1, const uint32_t n_cols1,
                                 const uint32_t n_cols2) {
    if (n_cols2 == 1) {
        mirath_matrix_ff_mu_mul_u256_vector(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }

    // set this to the number of element in each column*16 you want
    // to be able to load
    __m256i A[4];
    // number of avc register in each column in A
    const uint32_t nr256= n_rows1 >> 4u;
    const uint32_t tb   = n_rows1 - (nr256 * 16);
    const uint8_t tail  = (nr256 == 0) || (tb != 0);
    uint16_t tmp[16] __attribute__((aligned(32)));
    for (uint32_t j = 0; j < n_cols1; ++j) {
        const uint16_t *m1 = (uint16_t*)(matrix1 + j*n_rows1);
        for (uint32_t i = 0; i < nr256; i++) {
            A[i] = _mm256_loadu_si256((const __m256i *)(m1 + i*16));
        }
        if (tail) {
            for (uint32_t k = 0; k < tb; ++k) {
                tmp[k] = m1[nr256*16 + k];
            }
            A[nr256] = _mm256_load_si256((const __m256i *)tmp);
        }
        // iterate over one row in matrix2
        for (uint32_t i = 0; i < n_cols2; ++i) {
            const __m256i b = _mm256_set1_epi16(matrix2[i*n_cols1 + j]);
            for (uint32_t s = 0; s < nr256; ++s) {
                __m256i r = mirath_ff_mu_mul_u256(A[s], b);
                r ^= _mm256_loadu_si256((__m256i *)(result + i*n_rows1 + s*16));
                _mm256_storeu_si256((__m256i *)(result + i*n_rows1 + s*16), r);
            }
            // tail mngt
            if (tail) {
                const __m256i r = mirath_ff_mu_mul_u256(A[nr256], b);
                _mm256_store_si256((__m256i *)tmp, r);
                for (uint32_t k = 0; k < tb; ++k) {
                    result[i*n_rows1 + nr256*16 + k] ^= tmp[k];
                }
            }
        }
    }
}

/**
 * \fn static inline void mirath_matrix_ff_mu_add_product(ff_mu_t *result, const ff_mu_t *matrix1,
 *                                   const ff_mu_t *matrix2, const uint32_t n_rows1,
 *                                   const uint32_t n_cols1, const uint32_t n_cols2)
 * \brief result += matrix1 * matrix2
 *
 * \param[out] result Matrix over ff_mu
 * \param[in] matrix1 Matrix over ff_mu
 * \param[in] matrix2 Matrix over ff_mu
 * \param[in] n_rows1 number of rows in matrix1
 * \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
 * \param[in] n_cols2 number of columns in matrix2
 */
static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1,
                                     const ff_mu_t *matrix2, const uint32_t n_rows1,
                                     const uint32_t n_cols1, const uint32_t n_cols2) {
    memset(result, 0, n_rows1 * n_cols2 * sizeof(ff_mu_t));
    mirath_matrix_ff_mu_add_product(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
}

/**
 * \fn static inline void mirath_matrix_map_ff_to_ff_mu(ff_mu_t *out, const uint8_t *input, const uint32_t nrows, const uint32_t ncols)
 * \brief mapping from ff to ff_mu
 *
 * \param[out] out Matrix over ff_mu
 * \param[in] input Matrix over ff
 * \param[in] nrows number of rows
 * \param[in] ncols number of columns
 */
static inline void mirath_matrix_map_ff_to_ff_mu(ff_mu_t *out, const uint8_t *input, const uint32_t nrows, const uint32_t ncols) {
    if ((nrows % 8 == 0) || (ncols == 1)) {
        mirath_vector_ff_mu_set_to_ff(out, input, nrows*ncols);
        return;
    }

    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(nrows);
    for (uint32_t col = 0; col < ncols; ++col) {
        const ff_t *in2 = input + col*gf2_col_bytes;
        ff_mu_t *o = out + col*nrows;

        mirath_vector_ff_mu_set_to_ff(o, in2, nrows);
    }
}

#endif
