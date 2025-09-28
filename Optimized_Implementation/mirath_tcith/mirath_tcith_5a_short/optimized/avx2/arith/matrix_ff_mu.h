#ifndef MIRATH_MATRIX_FF_MU_H
#define MIRATH_MATRIX_FF_MU_H

#include <stdint.h>
#include <string.h>

#include "ff_mu.h"
#include "ff.h"
#include "matrix_ff_arith.h"
#include "vector_ff_mu.h"

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

/// \param m
/// \param nrows
/// \param i
/// \param j
/// \param bytes
/// \return
static inline
__m128i gf16to3_matrix_load8(const ff_mu_t *m,
                              const uint32_t nrows,
                              const uint32_t i,
                              const uint32_t j,
                              const uint32_t bytes) {
    uint8_t tmp[16] __attribute__((aligned(32)));

    const uint32_t pos = j*nrows + i;
    uint8_t *p = (uint8_t *)(m + pos);
    for (uint32_t k = 0; k < MIN(bytes, 16); k++) {
        tmp[k] = p[k];
    }

    return _mm_load_si128((__m128i *)tmp);
}


/// loads min(16, a)
/// \param m
/// \param a
/// \return
static inline
__m256i gf16to3_matrix_load16(const ff_mu_t *m,
                              const uint32_t a) {

    uint16_t tmp[16] __attribute__((aligned(32)));

    for (uint32_t k = 0; k < MIN(a, 16); k++) {
        tmp[k] = m[k];
    }

    return _mm256_load_si256((__m256i *)tmp);
}

///
/// \param m
/// \param t
/// \param a
static inline
void gf16to3_matrix_write16(ff_mu_t *m,
                            const __m256i t,
                            const uint32_t a) {
    uint16_t tmp[16] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i *)tmp, t);
    for (uint32_t k = 0; k < MIN(a, 16); k++) {
        m[k] = tmp[k];
    }
}

///
/// \param m
/// \param nrows
/// \param i
/// \param j
/// \param bytes
/// \param a
static inline
void gf16to3_matrix_store16(ff_mu_t *m,
                            const uint32_t nrows,
                            const uint32_t i,
                            const uint32_t j,
                            const uint32_t bytes,
                            const __m128i a) {
    uint8_t tmp[16] __attribute__((aligned(32)));
    _mm_store_si128((__m128i *)tmp, a);

    const uint32_t pos = j*nrows + i;
    uint8_t *p = (uint8_t *)(m + pos);
    for (uint32_t k = 0; k < MIN(bytes, 16); k++) {
        p[k] = tmp[k];
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
    if ((n_rows % 8 == 0) || (n_cols == 1)) {
        mirath_vector_ff_mu_add_ff(matrix1, matrix2, matrix3, n_rows*n_cols);
        return;
    }

    for (uint32_t j = 0; j < n_cols; j++) {
        uint32_t i = 0;
        const uint32_t off = j*((n_rows+1)/2);
        for (; (i+16) <= n_rows; i+= 16) {
            const uint32_t t11 = *((uint32_t *)(matrix3 + off + i + 0));
            const uint32_t t12 = *((uint32_t *)(matrix3 + off + i + 4));
            const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const uint64_t t22 = _pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);
            const __m128i t1 = _mm_set_epi64x(t22, t21);
            const __m256i m3 = _mm256_cvtepu8_epi16(t1);

            const __m256i m2 = _mm256_loadu_si256((const __m256i *)(matrix2 + j*n_rows + i*2));

            const __m256i m1 = m3 ^ m2;
            _mm256_storeu_si256((__m256i *)(matrix1 + j*n_rows + i*2), m1);
        }

        for (; (i+8) <= n_rows; i+= 8) {
            const uint32_t t11 = *((uint32_t *)(matrix3 + off + i + 0));
            const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const __m128i t1 = _mm_set_epi64x(t21, t21);
            const __m128i m3 = _mm_cvtepu8_epi16(t1);

            const __m128i m2 = _mm_loadu_si128((const __m128i *)(matrix2 + j*n_rows + i*2));

            const __m128i m1 = m3 ^ m2;
            _mm_storeu_si128((__m128i *)(matrix1 + j*n_rows + i*2), m1);
        }

        if (i < n_rows) {
            const uint32_t rbytes = (n_rows - i + 1) / 2;
            const uint32_t wbytes = (n_rows - i) * 2;
            const uint32_t t11 = gf16_matrix_load4(matrix3, n_rows, i, j, rbytes);
            const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const __m128i t1 = _mm_set_epi64x(t21, t21);
            const __m128i m3 = _mm_cvtepu8_epi16(t1);

            // const __m128i m2 = _mm_loadu_si128((const __m128i *)(matrix2 + j*n_rows + i*2));
            const __m128i m2 = gf16to3_matrix_load8(matrix2, n_rows, i, j, wbytes);
            const __m128i m1 = m3 ^ m2;
            gf16to3_matrix_store16(matrix1, n_rows, i, j, wbytes, m1);
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
    if ((n_rows%2) == 0) {
        mirath_vector_ff_mu_add_multiple_ff(matrix1, matrix1, scalar, matrix2, n_rows * n_cols);
        return;
    }

    // = {0} to satisfy valgrind
    uint8_t tmp[32] __attribute__((aligned(32))) = {0};
    uint16_t *tmp2 = (uint16_t *)tmp;

    const __m256i s = _mm256_set1_epi16(scalar);
    const uint32_t bytes = (n_rows + 1) / 2;
    const uint32_t limit = n_rows % 16;

    for (uint32_t k = 0; k < n_cols; k++) {
        uint32_t i = n_rows;

        ff_mu_t *r    = matrix1 + k*n_rows;
        const ff_t *m = matrix2 + k*bytes;

        // avx2 code
        while (i >= 16u) {
            const __m256i t0 = mirath_vector_extend_gf16_x16(m);
            const __m256i t1 = t0 ^ _mm256_slli_epi16(t0, 8);
            const __m256i t2 = mirath_ff_mul_u256(s, t1);
            const __m256i t3 = _mm256_loadu_si256((const __m256i *)r);

            _mm256_storeu_si256((__m256i *)r, t2^t3);
            i -= 16u;
            m += 8u;
            r += 16u;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { tmp[j] = m[j]; }

            const __m256i t0 = mirath_vector_extend_gf16_x16(tmp);
            const __m256i t1 = t0 ^ _mm256_slli_epi16(t0, 8);
            const __m256i t2 = mirath_ff_mul_u256(s, t1);

            _mm256_store_si256((__m256i *)tmp, t2);

            for (uint32_t j = 0; j < limit; ++j) { r[j] ^= tmp2[j]; }
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
/// Assumes that matrix2 is a vector
/// \param[out] result Matrix over gf16to3
/// \param[in] matrix1 Matrix over gf16
/// \param[in] matrix2 Matrix over gf16to3
/// \param[in] n_rows number of rows in matrix1
/// \param[in] n_cols number of columns and rows in matrix1 and matrix2 respectively
static inline void mirath_matrix_ff_mu_product_vector_ff1mu(ff_mu_t *result,
                                                          const ff_t *matrix1,
                                                          const ff_mu_t *matrix2,
                                                          const uint32_t n_rows,
                                                          const uint32_t n_cols) {
    const uint32_t bytes = (n_rows+1)/2;
    const uint32_t limit = n_rows % 16;

    for (uint32_t k = 0; k < n_cols; k++) {
        uint32_t i = n_rows;

        const __m256i b = _mm256_set1_epi16(matrix2[k] & 0XFFF);
        ff_mu_t    *r = result;
        const ff_t *m = matrix1 + k*bytes;

        while (i >= 16) {
            const uint32_t t11 = *((uint32_t *)(m + 0));
            const uint32_t t12 = *((uint32_t *)(m + 4));
            const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const uint64_t t22 = _pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);
            const __m128i t1 = _mm_set_epi64x(t22, t21);
            const __m256i a1 = _mm256_cvtepu8_epi16(t1);
            const __m256i a2 = a1 ^ _mm256_slli_epi16(a1, 8);
            const __m256i t2 = mirath_ff_mul_u256(b, a2);
            const __m256i r1 = _mm256_loadu_si256((const __m256i *)r);
            const __m256i t3 = t2 ^ r1;
            _mm256_storeu_si256((__m256i *)r, t3);

            m += 8;
            r += 16;
            i -= 16;
        }

        if (limit) {
            ff_t tmp[32] __attribute__((aligned(32))) = {0};
            uint32_t *tmp2 = (uint32_t *)tmp;
            uint16_t *tmp3 = (uint16_t *)tmp;
            for (uint32_t j = 0; j < ((limit + 1)/2); j++) { tmp[j] = m[j];}

            const uint32_t t11 = *(tmp2 + 0);
            const uint32_t t12 = *(tmp2 + 1);
            const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
            const uint64_t t22 = _pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);
            const __m128i t1 = _mm_set_epi64x(t22, t21);
            const __m256i a1 = _mm256_cvtepu8_epi16(t1);
            const __m256i a2 = a1 ^ _mm256_slli_epi16(a1, 8);
            const __m256i t2 = mirath_ff_mul_u256(b, a2);
            _mm256_store_si256((__m256i *)tmp, t2);

            for (uint32_t j = 0; j < limit; j++) { r[j] ^= tmp3[j];}
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
        memset(result, 0, n_rows1*n_cols2*sizeof(ff_mu_t));
        mirath_matrix_ff_mu_product_vector_ff1mu(result, matrix1, matrix2, n_rows1, n_cols1);
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

/// assumes that B is a column vector
/// \param result
/// \param A
/// \param B
/// \param nrows in A
/// \param ncols in A
static inline void gf16to3_matrix_mul_vector(ff_mu_t *result,
                                             const ff_mu_t *A,
                                             const ff_mu_t *B,
                                             const uint32_t nrows,
                                             const uint32_t ncols) {
    const uint32_t limit = nrows % 16;
    const __m256i mask = _mm256_set1_epi16(0xFFF);
    for (uint32_t col = 0; col < ncols; ++col) {
        ff_mu_t        *r = result;
        const ff_mu_t *m1 = A + col*nrows;
        const __m256i b1 = _mm256_set1_epi16(B[col]);
        uint32_t i = nrows;
        while (i >= 16) {
            const __m256i a1 = _mm256_loadu_si256((const __m256i *)m1) & mask;
            const __m256i c1 = mirath_ff_mu_mul_u256(a1, b1);

            const __m256i r1 = _mm256_loadu_si256((const __m256i *)r);
            const __m256i c2 = c1 ^ r1;
            _mm256_storeu_si256((__m256i *)r, c2);

            i  -= 16;
            r  += 16;
            m1 += 16;
        }

        if (limit) {
            const __m256i a1 = gf16to3_matrix_load16(m1, limit) & mask;
            const __m256i c1 = mirath_ff_mu_mul_u256(a1, b1);
            const __m256i r1 = gf16to3_matrix_load16(r, limit);
            const __m256i c2 = c1 ^ r1;
            gf16to3_matrix_write16(r, c2, limit);
        }
    }
}

/// Assumes 16 < rows <= 32 in A, the rest doesnt matter
/// col major
static inline void gf16to3_matrix_mul_le32xCxC(ff_mu_t*result,
                                               const ff_mu_t*A,
                                               const ff_mu_t*B,
                                               const uint16_t nrows,
                                               const uint32_t nrows_B,
                                               const uint32_t ncols_B) {
    const uint32_t tail = nrows % 16;
    for (uint32_t j = 0; j < nrows_B; ++j) {
        const __m256i a1 = _mm256_loadu_si256((const __m256i *)(A + j*nrows));
        const __m256i a2 = gf16to3_matrix_load16(A + j*nrows + 16, tail);

        for (uint32_t i = 0; i < ncols_B; ++i) {
            const __m256i b = _mm256_set1_epi16(*(B + i*nrows_B + j));
            const __m256i c1 = mirath_ff_mu_mul_u256(a1, b);
            const __m256i c2 = mirath_ff_mu_mul_u256(a2, b);

            const __m256i r1 = _mm256_loadu_si256((const __m256i *)(result + i*nrows));
            const __m256i r2 = gf16to3_matrix_load16(result + i*nrows + 16, tail);
            const __m256i t1 = c1 ^ r1;
            const __m256i t2 = c2 ^ r2;

            _mm256_storeu_si256((__m256i *)(result + i*nrows), t1);
            gf16to3_matrix_write16(result + i*nrows + 16, t2, tail);
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
    memset(result, 0, n_rows1*n_cols2*sizeof(ff_mu_t));
    if (n_cols2 == 1) {
        gf16to3_matrix_mul_vector(result, matrix1, matrix2, n_rows1, n_cols1);
    } else {
        gf16to3_matrix_mul_le32xCxC(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
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

