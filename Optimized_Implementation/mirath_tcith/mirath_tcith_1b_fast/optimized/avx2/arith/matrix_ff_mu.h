#ifndef MIRATH_MATRIX_FF_MU_8_H
#define MIRATH_MATRIX_FF_MU_8_H

#include <stdint.h>
#include <stdio.h>

#include "matrix_ff_arith.h"
#include "vector_ff_mu.h"
#include "ff_mu.h"

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


/// reads up to 16 bytes from a single columns
static inline __m128i read16column(const ff_mu_t*col,
                                   const uint32_t bytes) {
    if (bytes == 32) { return _mm_loadu_si128((const __m128i *)col); }

    uint8_t tmp[16] __attribute__((aligned(16)));
    for (uint32_t i = 0; i < bytes; ++i) {
        tmp[i] = col[i];
    }

    return _mm_load_si128((const __m128i *)tmp);
}

/// reads up to 32 bytes from a single columns
static inline __m256i read32column(const ff_mu_t*col,
                                   const uint32_t bytes) {
    if (bytes == 32) { return _mm256_loadu_si256((const __m256i *)col); }

    uint8_t tmp[32] __attribute__((aligned(32)));
    for (uint32_t i = 0; i < bytes; ++i) {
        tmp[i] = col[i];
    }

    return _mm256_load_si256((const __m256i *)tmp);
}

/// writes up to 32 bytes into a single colum
static inline void write32column(ff_mu_t*col,
                                 const __m256i data,
                                 const uint32_t bytes) {
    uint8_t tmp[32] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i *)tmp, data);
    for (uint32_t i = 0; i < bytes; ++i) {
        col[i] = tmp[i];
    }
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
    mirath_vector_ff_mu_add(matrix1, matrix2, matrix3, n_rows * n_cols);
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
    if ((n_rows % 8 == 0) || (n_cols == 1)) {
        mirath_vector_ff_mu_add_ff(matrix1, matrix2, matrix3, n_rows*n_cols);
        return;
    }

    const uint32_t bytes_per_column = mirath_matrix_ff_bytes_per_column(n_rows);
    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = n_rows;
        ff_t *in   = (ff_t *)matrix3 + j*bytes_per_column;
        ff_mu_t *in2 = (ff_mu_t *)matrix2 + j*n_rows;
        ff_mu_t *out = matrix1 + j*n_rows;
        while (i >= 32u) {
            const __m256i m1 = mirath_expand_gf2_x32_u256(in);
            const __m256i m2 = _mm256_loadu_si256((const __m256i *)in2);
            const __m256i t1 = m1 ^ m2;
            _mm256_storeu_si256((__m256i *)out, t1);
            in  += 4u;
            in2 += 32u;
            out += 32u;
            i   -= 32u;
        }

        while (i >= 16u) {
            const __m128i m1 = mirath_expand_gf2_x16_u256(in);
            const __m128i m2 = _mm_loadu_si128((const __m128i *)in2);
            const __m128i t1 = m1 ^ m2;
            _mm_storeu_si128((__m128i *)out, t1);
            in  += 2u;
            out += 16u;
            in2 += 16u;
            i   -= 16u;
        }

        while (i >= 8u) {
            const uint64_t m1 = mirath_expand_gf2_x8_u256(in);
            const uint64_t m2 = *(uint64_t *)(in2);
            const uint64_t t1 = m1 ^ m2;
            *(uint64_t *)out = t1;
            in  += 1u;
            out += 8u;
            in2 += 8u;
            i   -= 8u;
        }

        if (i) {
            uint8_t tmp[8] __attribute__((aligned(32)));
            const uint64_t m1 = mirath_expand_gf2_x8_u256(in);

            for (uint32_t k = 0; k < i; k++) { tmp[k] = in2[k]; }
            const uint64_t m2 = *(uint64_t *)(tmp);
            const uint64_t t1 = m1 ^ m2;
            *(uint64_t *)tmp = t1;
            for (uint32_t k = 0; k < i; k++) { out[k] = tmp[k]; }
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
 * \param[in] nrows number of rows
 * \param[in] ncols number of columns
 */
static inline void mirath_matrix_ff_mu_add_multiple_ff(ff_mu_t *matrix1, ff_mu_t scalar, const ff_t *matrix2,
                                         const uint32_t nrows, const uint32_t ncols) {
    if ((nrows % 8 == 0) || (ncols == 1)) {
        mirath_vector_ff_mu_add_multiple_ff(matrix1, matrix1, scalar, matrix2, nrows*ncols);
        return;
    }
    const __m128i s128 = _mm_set1_epi8(scalar);
    if (nrows == 4) {
        const uint64_t m = 0x01010101;
        const uint32_t limit = ncols % 4;

        const ff_t *mm2 = matrix2;
        ff_mu_t *out = matrix1;

        /// NOTE: asumes n_cols %2 == 2
        uint32_t col = 0;
        for (; (col+4) <= ncols; col+=4) {
            const uint32_t a = *(uint32_t *)(matrix2 + col);

            const uint64_t t41 = _pdep_u64(a>> 0u, m);
            const uint64_t t42 = _pdep_u64(a>> 8u, m);
            const uint64_t t43 = _pdep_u64(a>>16u, m);
            const uint64_t t44 = _pdep_u64(a>>24u, m);
            const uint64_t t31 =  t41 ^ (t42 << 32);
            const uint64_t t32 =  t43 ^ (t44 << 32);
            const __m128i t  = _mm_set_epi64x(t32, t31);

            const __m128i m1 = _mm_loadu_si128((__m128i *)out);
            const __m128i c1 = mirath_ff_mu_mul_ff_u128(s128, t);
            _mm_storeu_si128((__m128i *)(out +  0), m1 ^ c1);

            out += 4*nrows;
            mm2 += 4;
        }

        if (limit) {
            ff_mu_t tmp[16] __attribute__((aligned(32)));

            const uint64_t t41 = _pdep_u64(mm2[0], m);
            const uint64_t t42 = _pdep_u64(mm2[1], m);
            const uint64_t t31 = t41 ^ (t42 << 32);
            const __m128i t  = _mm_set_epi64x(0, t31);
            const __m128i c1 = mirath_ff_mu_mul_ff_u128(s128, t);
            _mm_store_si128((__m128i *)tmp, c1);

            for (uint32_t i = 0; i < 2*nrows; i++) {
                out[i] ^= tmp[i];
            }
        }
        return;
    }


    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(nrows);
    for (uint32_t col = 0; col < ncols; ++col) {
        ff_mu_t *o     = matrix1 + col*nrows;
        const ff_t *in = matrix2 + col*gf2_col_bytes;

        mirath_vector_ff_mu_add_multiple_ff(o, o, scalar, in, nrows);
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
/// `n_cols2` == 1
///
/// \param[out] result Matrix over gf256
/// \param[in] matrix1 Matrix over gf2
/// \param[in] matrix2 Matrix over gf256
/// \param[in] n_rows1 number of rows in matrix1
/// \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
static inline void mirath_matrix_ff_mu_product_ff1mu_vector(ff_mu_t *result,
                                                            const ff_t *matrix1,
                                                            const ff_mu_t *matrix2,
                                                            const uint32_t n_rows1,
                                                            const uint32_t n_cols1) {
    uint8_t tmp[32] __attribute__((aligned(32)));

    const uint32_t limit = n_rows1 % 32;
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col * bytes_per_col;
        const __m256i b = _mm256_set1_epi8(*(matrix2 + col));
        ff_mu_t *r = result;
        while ((i + 32) <= n_rows1) {
            const __m256i a = mirath_expand_gf2_x32_u256(m1);
            __m256i t = mirath_ff_mu_mul_ff_u256(b, a);
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 4;
            r  += 32;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+7)/8; ++j) { tmp[j] = m1[j]; }
            const __m256i a = mirath_expand_gf2_x32_u256(tmp);
            __m256i t = mirath_ff_mu_mul_ff_u256(b, a);
            _mm256_store_si256((__m256i *)tmp, t);
            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp[j]; }
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
        memset(result, 0, n_rows1);
        mirath_matrix_ff_mu_product_ff1mu_vector(result, matrix1, matrix2, n_rows1, n_cols1);
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

/// special matrix multiplication
/// B <= 8
/// 32 < a <= 64
/// \param[out] C Matrix over ff_mu_t a by c matrix
/// \param[in]  A Matrix over ff_mu_t a by b matrix
/// \param[in]  B Matrix over ff_mu_t b by c matrix
static inline void mirath_matrix_product_le64xBxle48_u256(ff_mu_t*C,
                                                         const ff_mu_t*A,
                                                         const ff_mu_t*B,
                                                         const uint32_t a,
                                                         const uint32_t b,
                                                         const uint32_t c) {
    const __m256i m = _mm256_set1_epi8(0x0F);
    const uint32_t la = a % 32;

    // for each col in A
    for (uint32_t i = 0; i < b; ++i) {
        // load full col in A
        const __m256i col_A1 = _mm256_loadu_si256((const __m256i *)(A + i*a));
        const __m256i col_A2 = read32column(A + i*a + 32, la);

        // for each element in row i in B
        for (uint32_t j = 0; j < c; ++j) {
            // load element in B
            const ff_mu_t elm_B = B[j*b + i];
            const __m256i m_tab = mirath_ff_mu_generate_multab_16_single_element_u256(elm_B);
            const __m256i ml = _mm256_permute2x128_si256(m_tab, m_tab, 0);
            const __m256i mh = _mm256_permute2x128_si256(m_tab, m_tab, 0x11);

            const __m256i r1 = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, col_A1, m);
            const __m256i r2 = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, col_A2, m);


            const __m256i c1 = _mm256_loadu_si256((const __m256i *)(C + j*a));
            const __m256i c2 = read32column(C + j*a + 32, la);
            _mm256_storeu_si256((__m256i *)(C + j*a), c1^r1);
            write32column(C + j*a + 32, c2^r2, la);
        }
    }
}



/// \brief result = matrix1 * matrix2
/// `n_cols2` == 1
///
/// \param[out] result Matrix over gf256
/// \param[in] matrix1 Matrix over gf2
/// \param[in] matrix2 Matrix over gf256
/// \param[in] n_rows1 number of rows in matrix1
/// \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
static inline void mirath_matrix_product_u256_vector(ff_mu_t*result,
                                                     const ff_mu_t*matrix1,
                                                     const ff_mu_t*matrix2,
                                                     const uint32_t n_rows1,
                                                     const uint32_t n_cols1) {
    uint8_t tmp[32] __attribute__((aligned(32)));

    const uint32_t limit = n_rows1 % 32;

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const ff_mu_t*m1 = matrix1 + col*n_rows1;
        const __m256i b = _mm256_set1_epi8(*(matrix2 + col));
        ff_mu_t*r = result;
        while ((i + 32) <= n_rows1) {
            const __m256i a = _mm256_loadu_si256((const __m256i *)m1);
            __m256i t = gf256v_mul_u256(a, b);
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 32;
            r  += 32;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < limit; ++j) { tmp[j] = m1[j]; }
            const __m256i a = _mm256_loadu_si256((const __m256i *)tmp);
            __m256i t = gf256v_mul_u256(a, b);
            _mm256_store_si256((__m256i *)tmp, t);
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
static inline void mirath_matrix_ff_mu_product(ff_mu_t *result, const ff_mu_t *matrix1, const ff_mu_t *matrix2,
                                 const uint32_t n_rows1, const uint32_t n_cols1,
                                 const uint32_t n_cols2) {
    memset(result,0, n_rows1*n_cols2);
     if (n_cols2 == 1) {
        mirath_matrix_product_u256_vector(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }

    mirath_matrix_product_le64xBxle48_u256(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
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
    if (nrows % 8 == 0) {
        mirath_vector_ff_mu_set_to_ff_u256(out, input, nrows*ncols);
        return;
    }
    const uint32_t gf2_col_bytes = mirath_matrix_ff_bytes_per_column(nrows);
    for (uint32_t col = 0; col < ncols; ++col) {
        const ff_t *in2 = input + col*gf2_col_bytes;
        ff_mu_t *o = out + col*nrows;

        mirath_vector_ff_mu_set_to_ff_u256(o, in2, nrows);
    }
}

#endif
