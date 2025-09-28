#ifndef MATRIX_FF_ARITH_H
#define MATRIX_FF_ARITH_H

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

/// reads `bytes`from m + nrows*j/2 + i/2
/// \param m base pointer of the matrix
/// \param nrows number of rows
/// \param i number of column to read from
/// \param j number of row to read from
/// \param bytes number of bytes to read
static inline
uint32_t gf16_matrix_load4(const ff_t *m,
                           const uint32_t nrows,
                           const uint32_t i,
                           const uint32_t j,
                           const uint32_t bytes) {
    uint8_t tmp[4];
    const uint32_t pos = mirath_matrix_ff_bytes_per_column(nrows) * j + i/2;
    for (uint32_t k = 0; k < MIN(bytes, 4); k++) {
        tmp[k] = m[pos + k];
    }

    return *((uint32_t*)tmp);
}

/// clears the last byte in each column if the number of rows is uneven
static inline void mirath_matrix_set_to_ff(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    if (n_rows & 1) {
        const uint32_t matrix_height =  mirath_matrix_ff_bytes_per_column(n_rows);
        const uint32_t matrix_height_x = matrix_height -  1;

        for (uint32_t i = 0; i < n_cols; i++) {
            matrix[i * matrix_height + matrix_height_x ] &= 0x0f;
        }
    }
}

/// \return m[i, j]
static inline ff_t mirath_matrix_ff_get_entry(const ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    if (i & 1u) { // i is odd
        return  matrix[nbytes_col * j + (i >> 1)] >> 4;
    }
    else {
        return matrix[nbytes_col * j + (i >> 1)] & 0x0f;
    }
}

///  m[i, j] = scalar
static inline void mirath_matrix_ff_set_entry(ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j, const ff_t scalar) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    const uint32_t target_byte_id = nbytes_col * j + (i >> 1);
    if (i & 1) // i is odd
    {
        matrix[target_byte_id] &= 0x0f;
        matrix[target_byte_id] |= (scalar << 4);
    }
    else {
        matrix[target_byte_id] &= 0xf0;
        matrix[target_byte_id] |= scalar;
    }
}

/// matrix1 = matrix2 + matrix3
static inline void mirath_matrix_ff_add_arith(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
		const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);

    for (uint32_t i = 0; i < n_bytes; i++) {
        matrix1[i] = matrix2[i] ^ matrix3[i];
    }
}

/// matrix1 += scalar *matrix2
static inline void mirath_matrix_ff_add_multiple_arith(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);

    for (uint32_t i = 0; i < n_bytes; i++) {
        matrix1[i] ^= mirath_ff_product(scalar, matrix2[i] & 0xf);
        matrix1[i] ^= mirath_ff_product(scalar, matrix2[i] >> 4) << 4;
    }
}

/// \param result
/// \param matrix1
/// \param matrix2
/// \param n_rows1
/// \param n_cols1
static inline void mirath_matrix_ff_product_vector_u256(ff_t *result,
                                                        const ff_t *matrix1,
                                                        const ff_t *matrix2,
                                                        const uint32_t n_rows1,
                                                        const uint32_t n_cols1) {
    uint8_t tmp[32] __attribute__((aligned(32))) = {0};

    const uint32_t limit = n_rows1 % 64;
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);
    const __m256i mask = _mm256_set1_epi8(0xf);

    for (uint32_t col = 0; col < n_cols1; ++col) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col * bytes_per_col;
        const uint8_t b2 = mirath_matrix_ff_get_entry(matrix2, 0, col, 0);
        const __m256i bl = mirath_ff_tbl32_multab(b2);
        const __m256i bh = _mm256_slli_epi16(bl, 4);
        ff_t *r = result;

        while ((i + 64) <= n_rows1) {
            const __m256i a = _mm256_loadu_si256((__m256i *)m1);
            __m256i t = mirath_ff_linear_transform_8x8_256b(bl, bh, a, mask);

            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 32;
            r  += 32;
            i  += 64;
        }

        if (limit) {
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { tmp[j] = m1[j]; }
            const __m256i a = _mm256_loadu_si256((__m256i *)tmp);
            __m256i t = mirath_ff_linear_transform_8x8_256b(bl, bh, a, mask);
            _mm256_store_si256((__m256i *)tmp, t);
            for (uint32_t j = 0; j < (limit+1)/2; ++j) { r[j] ^= tmp[j]; }
        }
    }
}

/// \param c[out]: output matrix of size 16xk
/// \param a[in]: input matrix of size 16x4 
/// \param b[in]: input matrix of size 4xk
/// \param k[in]: number of columns in b
static inline
void mirath_matrix_ff_mul_16x4xk(uint8_t *__restrict c,
                                 const uint8_t *__restrict__ a,
                                 const uint8_t *__restrict__ b, 
                                 const uint32_t k) {
    // fully load A into a register 
    const __m256i A = _mm256_loadu_si256((const __m256i *)a);
    uint64_t *c64 = (uint64_t *)c;
    const uint32_t bytes_per_col_B = 2;

    // some helper variables
    const __m256i one = _mm256_set1_epi8(0x11);
    for (uint32_t i = 0; i < k; i++) {
        // extend a single column to an avx register
        const uint16_t t1 = *((uint16_t *)(b + i*bytes_per_col_B));
        const uint64_t t2 = _pdep_u64(t1, 0xF000F000F000F);
        const __m128i t3 = _mm_set_epi64x(0, t2);
        const __m256i t4 = _mm256_cvtepu16_epi64(t3);
        const __m256i t5 = _mm256_mul_epu32(t4, one);
        const __m256i t6 = t5^_mm256_slli_si256(t5, 4);

        const __m256i c1 = mirath_ff_mul_full_u256(t6, A);
        const uint64_t c2 = mirath_hadd_u256_64(c1);

        c64[i] = c2;
    }
}

/// result = matrix1 * matrix2
// matrix1 of size n_rows1 * n_cols1
// matrix2 of size n_cols1 * n_cols2
// result  of size n_rows1 * n_cols2
static inline void mirath_matrix_ff_product_arith(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    if ((n_rows1 == 16) && (n_cols1 == 4)) {
        mirath_matrix_ff_mul_16x4xk(result, matrix1, matrix2, n_cols2);
        return;
    }

    if (n_cols2 == 1) {
        mirath_matrix_ff_product_vector_u256(result, matrix1, matrix2, n_rows1, n_cols1);
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
