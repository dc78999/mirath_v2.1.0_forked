#ifndef MATRIX_FF_ARITH_H
#define MATRIX_FF_ARITH_H

#include <string.h>
#include "ff.h"
#include "vector_ff_arith.h"

#define mirath_matrix_ff_bytes_per_column(n_rows) (((n_rows) >> 3) + ((n_rows % 8) == 0 ? 0 : 1))
#define mirath_matrix_ff_bytes_size(n_rows, n_cols) ((mirath_matrix_ff_bytes_per_column(n_rows)) * (n_cols))

#define MIRATH_VAR_FF_AUX_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R) + mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_S_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_R))
#define MIRATH_VAR_FF_C_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_H_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K))
#define MIRATH_VAR_FF_Y_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, 1))
#define MIRATH_VAR_FF_T_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M, MIRATH_PARAM_N - MIRATH_PARAM_R))
#define MIRATH_VAR_FF_E_BYTES (mirath_matrix_ff_bytes_size(MIRATH_PARAM_M * MIRATH_PARAM_N, 1))

#define OFF_E_A ((8 * MIRATH_VAR_FF_Y_BYTES) - (MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K))
#define OFF_E_B ((8 * mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1)) - MIRATH_PARAM_K)

/// read `len` bytes from `ptr`
/// \param ptr[in]: pointer to read from
/// \param len[in]: number of bytes to read
static inline uint64_t gf2_read8rows(const uint8_t *ptr,
                                     const uint32_t len) {
    uint64_t t = 0;
    for (uint32_t i = 0; i < len; i++) {
        t ^= ((uint64_t )ptr[i]) << (i*8);
    }

    return t;
}

/// clears the last byte in each column, if its not completly used
/// \param matrix[in/out]: the matrix to clear
/// \param n_rows[in]: number of rows in matrix
/// \param n_cols[int]: number of rows i matrix
static inline void mirath_matrix_set_to_ff(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    if (n_rows & 0x7) {
        const uint32_t matrix_height =  mirath_matrix_ff_bytes_per_column(n_rows);
        const uint32_t matrix_height_x = matrix_height -  1;

        ff_t mask = 0xff >> (8 - (n_rows % 8));

        for (uint32_t i = 0; i < n_cols; i++) {
            matrix[i * matrix_height + matrix_height_x ] &= mask;
        }
    }
}

/// \param matrix[in]: matrix to access
/// \param n_rows[in]: number of rows in the matrix
/// \param i row[in]; to set
/// \param j column[in]: to set
/// \return m[i, j]
static inline ff_t mirath_matrix_ff_get_entry(const ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    const uint32_t idx_line = i / 8;
    const uint32_t bit_line = i % 8;

    return (matrix[nbytes_col * j + idx_line] >> bit_line) & 0x01;
}

/// m[i, j] = scalar
/// \param matrix[in/out]: matrix to alter
/// \param n_rows[in]: number of rows in matrix
/// \param i[in]: row 
/// \param j[in]: column
/// \param scalar[in]: input value
static inline void mirath_matrix_ff_set_entry(ff_t *matrix, const uint32_t n_rows, const uint32_t i, const uint32_t j, const ff_t scalar) {
    const uint32_t nbytes_col = mirath_matrix_ff_bytes_per_column(n_rows);
    const uint32_t idx_line = i / 8;
    const uint32_t bit_line = i % 8;

    const uint8_t mask = 0xff ^ (1 << bit_line);

    matrix[nbytes_col * j + idx_line] = (matrix[nbytes_col * j + idx_line] & mask) ^ (scalar << bit_line);
}

/// matrix1 = matrix2 + matrix3
static inline void mirath_matrix_ff_add_arith(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
		const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    mirath_vec_ff_add_arith(matrix1, matrix2, matrix3, n_bytes);
}

/// matrix1 += scalar *matrix2
static inline void mirath_matrix_ff_add_multiple_arith(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    mirath_vec_ff_add_multiple_arith(matrix1, matrix1, scalar, matrix2, n_bytes);
}


/// NOTE: assumes that `n_cols2` = 1
/// result = matrix1 * matrix2
/// matrix1 of size n_rows1 * n_cols1
/// matrix2 of size n_cols1 * n_cols2=1
/// result  of size n_rows1 * n_cols2=1
static inline void mirath_matrix_ff_mul_u256_vector(ff_t *result, 
                                                    const ff_t *matrix1, 
                                                    const ff_t *matrix2,
                                                    const uint32_t n_rows1, 
                                                    const uint32_t n_cols1) {
    uint8_t tmp[32] __attribute__((aligned(32)));

    const uint32_t limit = n_rows1 % 256;
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);

    for (uint32_t col = 0; col < n_cols1; col++) {
        uint32_t i = 0;
        const uint8_t *m1 = matrix1 + col*bytes_per_col;
        const char b_ = mirath_vec_ff_get_entry(matrix2, col);
        const __m256i b = _mm256_set1_epi8(-b_);
        uint8_t *r = result;

        while ((i + 256) <= n_rows1) {
            const __m256i a = _mm256_loadu_si256((const __m256i *)m1);
            __m256i t = a & b;
            t ^= _mm256_loadu_si256((__m256i *)(r));
            _mm256_storeu_si256((__m256i *)(r), t);

            m1 += 32;
            r  += 32;
            i  += 256;
        }
        
        if (limit) {
            for (uint32_t j = 0; j < (limit+7)/8; ++j) { tmp[j] = m1[j]; }
            const __m256i a = _mm256_loadu_si256((const __m256i *)tmp);
            __m256i t = a & b;
            _mm256_store_si256((__m256i *)tmp, t);
            for (uint32_t j = 0; j < (limit+7)/8; ++j) { r[j] ^= tmp[j]; }
        }
    }
}


/// NOTE: assumes n_cols1 <= 8
/// result = matrix1 * matrix2
/// matrix1 of size n_rows1 * n_cols1
/// matrix2 of size n_cols1 * n_cols2
/// result  of size n_rows1 * n_cols2
static inline void mirath_matrix_ff_mul_u256_Axle8xC(ff_t *result, 
                                                     const ff_t *matrix1, 
                                                     const ff_t *matrix2,
                                                     const uint32_t n_rows1, 
                                                     const uint32_t n_cols1, 
                                                     const uint32_t n_cols2) {
    const uint32_t bytes_per_col = mirath_matrix_ff_bytes_per_column(n_rows1);
    for (uint32_t colA = 0; colA < n_cols1; colA++) {
        const uint8_t *m1 = matrix1 + (colA*bytes_per_col);
        // load up to 8 columns
        uint64_t a = gf2_read8rows(m1, bytes_per_col);
        // iterate over all columns in B
        for (uint32_t col = 0; col < n_cols2; col++) {
            const uint64_t b1 = (matrix2[col]>>colA) & 1u;
            const uint64_t b = -b1;
            uint64_t c = a&b;
            
            uint8_t *r = result + col*bytes_per_col;
            for (uint32_t i = 0; i < bytes_per_col; i++) {
                r[i] ^= c;
                c >>= 8u;
            }
        }
    }
}

/// result = matrix1 * matrix2
/// matrix1 of size n_rows1 * n_cols1
/// matrix2 of size n_cols1 * n_cols2
/// result  of size n_rows1 * n_cols2
static inline void mirath_matrix_ff_product_arith(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    if (n_cols2 == 1) {
        memset(result, 0, mirath_matrix_ff_bytes_size(n_rows1, n_cols2));
        mirath_matrix_ff_mul_u256_vector(result, matrix1, matrix2, n_rows1, n_cols1);
        return;
    }

    if (n_cols1 <= 8) {
        memset(result, 0, mirath_matrix_ff_bytes_size(n_rows1, n_cols2));
        mirath_matrix_ff_mul_u256_Axle8xC(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
        return;
    }

    // fallback implementation

    uint32_t i, j, k;
    ff_t entry_i_k, entry_k_j, entry_i_j;

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
}

#endif
