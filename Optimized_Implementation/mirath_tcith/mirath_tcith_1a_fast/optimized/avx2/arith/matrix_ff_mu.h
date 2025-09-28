#ifndef MIRATH_MATRIX_FF_MU_H
#define MIRATH_MATRIX_FF_MU_H

#include <stdint.h>
#include <string.h>

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

/// matrix1 = matrix2
static inline void mirath_matrix_ff_mu_copy(ff_mu_t *matrix1, const ff_mu_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    memcpy(matrix1, matrix2, mirath_matrix_ff_mu_bytes_size(n_rows, n_cols));
}

/// reads up to 32 bytes from a single columns
/// \param col[in]: the column from which data is read
/// \param bytes[in]: number of bytes to read
static inline __m256i read32column(const ff_mu_t *col,
                                   const uint32_t bytes) {
    if (bytes == 32) { return _mm256_loadu_si256((const __m256i *)col); }

    uint8_t tmp[32] __attribute__((aligned(32)));
    for (uint32_t i = 0; i < bytes; ++i) {
        tmp[i] = col[i];
    }

    return _mm256_load_si256((const __m256i *)tmp);
}

/// writes up to "bytes" (<=32 bytes) into a single colum
/// \param col[in/out] `bytes` bytes will we written into memory
/// \param data[in]: the avx register to be written into memory
/// \param bytes[in]: number of bytes written to memory
static inline void write32column(ff_mu_t *col,
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
        mirath_vector_ff_mu_add_ff(matrix1, matrix2, matrix3, n_rows*n_cols);
        return;
    }

    const __m128i perm128 = _mm_load_si128((const __m128i *)mirath_map_ff_to_ff_mu);
    const uint32_t limit = n_rows % 16;

    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = n_rows;
        ff_t *in     = (ff_t    *)matrix3 + j*((n_rows + 1)/2);
        ff_mu_t *in2 = (ff_mu_t *)matrix2 + j*n_rows;
        ff_mu_t *out = matrix1 + j*n_rows;

        while (i >= 16u) {
            const __m128i m2 = _mm_loadu_si128((const __m128i *)in2);
            const uint32_t t11 = *(uint32_t *)(in +  0);
            const uint32_t t12 = *(uint32_t *)(in +  4);

            const uint64_t t21 = _pdep_u64(t11, 0x0f0f0f0f0f0f0f0f);
            const uint64_t t22 = _pdep_u64(t12, 0x0f0f0f0f0f0f0f0f);

            const __m128i t3 = _mm_setr_epi64((__m64)t21, (__m64)t22);
            const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
            const __m128i t5 = t4 ^ m2;
            _mm_storeu_si128((__m128i *)out, t5);

            in  += 8u;
            out += 16u;
            in2 += 16u;
            i   -= 16u;
        }

        if (limit) {
            uint8_t tmp1[16] __attribute__((aligned(32))) = {0};
            uint8_t tmp2[16] __attribute__((aligned(32))) = {0};
            for (uint32_t k = 0; k < i; k++) { tmp1[k] = in2[k]; }
            const __m128i m2 = _mm_load_si128((const __m128i *) tmp1);

            for (uint32_t k = 0; k < (i + 1) / 2; k++) { tmp2[k] = in[k]; }
            const __m64 t21 = (__m64) _pdep_u64(*(uint32_t *) (tmp2 + 0), 0x0f0f0f0f0f0f0f0f);
            const __m64 t22 = (__m64) _pdep_u64(*(uint32_t *) (tmp2 + 4), 0x0f0f0f0f0f0f0f0f);

            const __m128i t3 = _mm_setr_epi64(t21, t22);
            const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
            const __m128i t5 = t4 ^ m2;
            _mm_store_si128((__m128i *) tmp1, t5);
            for (uint32_t k = 0; k < i; k++) { out[k] = tmp1[k]; }
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
    if ((n_rows % 2) == 2) {
        mirath_vector_ff_mu_add_multiple_ff(matrix1, matrix1, scalar, matrix2, n_rows * n_cols);
        return;
    }

    const __m128i perm128 = _mm_load_si128((const __m128i *)mirath_map_ff_to_ff_mu);
    const __m256i tab = mirath_ff_mu_generate_multab_16_single_element_u256(scalar);
    const __m256i ml = _mm256_permute2x128_si256(tab, tab, 0);
    const __m256i mh = _mm256_permute2x128_si256(tab, tab, 0x11);

    const __m128i ml128 = _mm256_extracti128_si256(ml, 0);
    const __m128i mh128 = _mm256_extracti128_si256(mh, 0);
    const __m128i mask128 = _mm_set1_epi8(0xF);

    const uint32_t limit = n_rows % 16;

    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = n_rows;
        ff_t *in     = (ff_t    *)matrix2 + j*((n_rows + 1)/2);
        ff_mu_t *out = matrix1 + j*n_rows;

        while (i >= 16u) {
            const __m128i m2 = _mm_loadu_si128((const __m128i *)out);
            const uint32_t t11 = *(uint32_t *)(in +  0);
            const uint32_t t12 = *(uint32_t *)(in +  4);

            const __m64 t21 = (__m64)_pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const __m64 t22 = (__m64)_pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);

            const __m128i t3 = _mm_setr_epi64(t21, t22);
            const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
            const __m128i t5 = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, t4, mask128);
            const __m128i t6 = t5 ^ m2;
            _mm_storeu_si128((__m128i *)out, t6);
            in  += 8u;
            out += 16u;
            i   -= 16u;
        }

        if (limit) {
            uint8_t tmp1[16] __attribute__((aligned(32))) = {0};
            uint8_t tmp2[16] __attribute__((aligned(32))) = {0};
            for (uint32_t k = 0; k < limit; k++) { tmp1[k] = out[k]; }
            const __m128i m2 = _mm_load_si128((const __m128i *) tmp1);

            for (uint32_t k = 0; k < (limit + 1) / 2; k++) { tmp2[k] = in[k]; }
            const __m64 t21 = (__m64) _pdep_u64(*(uint32_t *) (tmp2 + 0), 0x0F0F0F0F0F0F0F0F);
            const __m64 t22 = (__m64) _pdep_u64(*(uint32_t *) (tmp2 + 4), 0x0F0F0F0F0F0F0F0F);

            const __m128i t3 = _mm_setr_epi64(t21, t22);
            const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
            const __m128i t5 = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, t4, mask128);
            const __m128i t6 = t5 ^ m2;
            _mm_store_si128((__m128i *) tmp1, t6);
            for (uint32_t k = 0; k < limit; k++) { out[k] = tmp1[k]; }
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
    for (uint32_t i = 0; i < n_rows; i++) {
        for (uint32_t j = 0; j < n_cols; j++) {
            const ff_mu_t entry1 = mirath_matrix_ff_mu_get_entry(matrix2, n_rows, i, j);
            const ff_mu_t entry2 = mirath_matrix_ff_mu_get_entry(matrix3, n_rows, i, j);
            const ff_mu_t entry3 = entry1 ^ mirath_ff_mu_mult(scalar, entry2);

            mirath_matrix_ff_mu_set_entry(matrix1, n_rows, i, j, entry3);
        }
    }
}
/// \brief result = matrix1 * matrix2
/// `n_cols2` == 1
///
/// \param[out] result Matrix over gf256
/// \param[in] matrix1 Matrix over gf16
/// \param[in] matrix2 Matrix over gf256
/// \param[in] n_rows1 number of rows in matrix1
/// \param[in] n_cols1 number of columns and rows in matrix1 and matrix2 respectively
static inline void mirath_matrix_product_gf16_1_vector_u256(ff_mu_t *result,
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
            const __m256i a = mirath_ff_mu_extend_gf16_x32(m1);
            __m256i t = gf256v_mul_u256(b, a);
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 16;
            r  += 32;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { tmp[j] = m1[j]; }
            const __m256i a = mirath_ff_mu_extend_gf16_x32(tmp);
            __m256i t = gf256v_mul_u256(b, a);
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
        mirath_matrix_product_gf16_1_vector_u256(result, matrix1, matrix2, n_rows1, n_cols1);
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
                entry_k_j = mirath_matrix_ff_get_entry(matrix2, n_cols1, k, j);
                entry_k_j = mirath_map_ff_to_ff_mu[entry_k_j];
                entry_i_j ^= mirath_ff_mu_mult(entry_i_k, entry_k_j);
            }

            mirath_matrix_ff_mu_set_entry(result, n_rows1, i, j, entry_i_j);
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
static inline void mirath_matrix_product_vector_u256(ff_mu_t *result,
                                                  const ff_mu_t *matrix1,
                                                  const ff_mu_t *matrix2,
                                                  const uint32_t n_rows1,
                                                  const uint32_t n_cols1) {
    uint8_t tmp[32] __attribute__((aligned(32)));

    const uint32_t limit = n_rows1 % 32;
    const __m256i m = _mm256_set1_epi8(0x0F);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const ff_mu_t *m1 = matrix1 + col*n_rows1;
        // const __m256i b = _mm256_set1_epi8(*(matrix2 + col));

        const __m256i m_tab = mirath_ff_mu_generate_multab_16_single_element_u256(matrix2[col]);
        const __m256i ml = _mm256_permute2x128_si256(m_tab, m_tab, 0);
        const __m256i mh = _mm256_permute2x128_si256(m_tab, m_tab, 0x11);

        ff_mu_t *r = result;
        while ((i + 32) <= n_rows1) {
            const __m256i a = _mm256_loadu_si256((const __m256i *)m1);
            // __m256i t = gf256v_mul_u256(a, b);
            __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, m);
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 32;
            r  += 32;
            i  += 32;
        }

        if (limit) {
            for (uint32_t j = 0; j < limit; ++j) { tmp[j] = m1[j]; }
            const __m256i a = _mm256_loadu_si256((const __m256i *)tmp);
            __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, m);
            _mm256_store_si256((__m256i *)tmp, t);
            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp[j]; }
        }
    }
}


/// special matrix multiplication
/// B <= 16
/// C <= 16
/// \param[out] C Matrix over ff_mu_t a by c matrix
/// \param[in]  A Matrix over ff_mu_t a by b matrix
/// \param[in]  B Matrix over ff_mu_t b by c matrix
/// \param[in] n_colsB number of columns in matrix2
static inline void mirath_matrix_product_le32xBxle16_u256(ff_mu_t *C,
                                                         const ff_mu_t *A,
                                                         const ff_mu_t *B,
                                                         const uint32_t a,
                                                         const uint32_t b,
                                                         const uint32_t c) {
    const __m256i m = _mm256_set1_epi8(0x0F);
    __m256i tmp_C[16] = {0};

    // for each col in A
    for (uint32_t i = 0; i < b; ++i) {
        // load full col in A
        const __m256i col_A = read32column(A + i*a, a);

        // for each element in row i in B
        for (uint32_t j = 0; j < c; ++j) {
            // load element in B
            const ff_mu_t elm_B = B[j*b + i];
            const __m256i m_tab = mirath_ff_mu_generate_multab_16_single_element_u256(elm_B);
            const __m256i ml = _mm256_permute2x128_si256(m_tab, m_tab, 0);
            const __m256i mh = _mm256_permute2x128_si256(m_tab, m_tab, 0x11);

            __m256i r = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, col_A, m);
            tmp_C[j] ^= r;
        }
    }

    for (uint32_t i = 0; i < c; ++i) {
        write32column(C +i*a, tmp_C[i], a);
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
    memset(result, 0, n_rows1*n_cols2);
    if (n_cols2 == 1) {
        mirath_matrix_product_vector_u256(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    } else {
        mirath_matrix_product_le32xBxle16_u256(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
        return;
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

