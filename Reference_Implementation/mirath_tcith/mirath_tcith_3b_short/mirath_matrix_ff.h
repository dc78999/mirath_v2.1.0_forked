#ifndef MIRATH_MATRIX_FF_H
#define MIRATH_MATRIX_FF_H

#include <stdint.h>
#include "hash_sha3.h"
#include "framework_data_types.h"
#include "mirath_arith.h"
#include "mirath_parameters.h"

/* Return the number of bytes of a 'n_rows x n_cols' matrix. */
//int matrix_bytes_size(int n_rows, int n_cols);

/* Initialized 'matrix' with zero entries. */
void mirath_matrix_init_zero(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols);

/* Initialized 'matrix' with random entries. */
void mirath_matrix_ff_init_random(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols, shake_prng_t *prng);

/* Generate the matrix 'H^prime' from the public seed 'seed_pk' */
void mirath_matrix_expand_seed_public_matrix(ff_t H[MIRATH_VAR_FF_H_BYTES], const seed_t seed_pk);

// /* Generate the matrices 'S' and 'C^prime' from the secret seed 'seed_sk' */
void mirath_matrix_expand_seed_secret_matrix(ff_t S[MIRATH_VAR_FF_S_BYTES], ff_t C[MIRATH_VAR_FF_C_BYTES], const seed_t seed_sk);

// /* Compute vector 'y' from 'H^prime', 'S', and 'C^prime'*/
void mirath_matrix_compute_y(ff_t y[MIRATH_VAR_FF_Y_BYTES],
                                    const ff_t S[MIRATH_VAR_FF_S_BYTES],
                                    const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                    const ff_t H[MIRATH_VAR_FF_H_BYTES]);

void mirath_tciht_compute_public_key(uint8_t *pk, const uint8_t *sk,
                                     const ff_t S[MIRATH_VAR_FF_S_BYTES], const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                     const ff_t H[MIRATH_VAR_FF_H_BYTES]);


void mirath_matrix_decompress_secret_key(ff_t S[MIRATH_VAR_FF_S_BYTES],
                      ff_t C[MIRATH_VAR_FF_C_BYTES],
                      ff_t H[MIRATH_VAR_FF_H_BYTES],
                      uint8_t *pk,
                      const uint8_t *sk);

void mirath_matrix_decompress_pk(ff_t H[MIRATH_VAR_FF_H_BYTES], ff_t y[MIRATH_VAR_FF_Y_BYTES], const uint8_t *pk);

/* Overwrite 'matrix1' with 'matrix2'. */
void mirath_matrix_ff_copy(ff_t *matrix1, const ff_t *matrix2, const uint32_t n_rows, const uint32_t n_cols);

/* Replace 'matrix' with '-matrix'. */
void mirath_matrix_ff_neg(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols);

/* set 'matrix1' with 'matrix2 + matrix3'. */
void mirath_matrix_ff_add(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 + scalar * matrix2'. */
void mirath_matrix_ff_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
                                   const uint32_t n_rows, const uint32_t n_cols);

/* set 'matrix1' with 'matrix2 - matrix3'. */
void mirath_matrix_ff_sub(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 - scalar * matrix2'. */
void mirath_matrix_ff_sub_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
                                   const uint32_t n_rows, const uint32_t n_cols);

/* Write 'matrix1 * matrix2' over 'result'. */
void mirath_matrix_ff_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
                              const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2);

/* Write '[matrix1 | matrix2]' over 'result'. */
void mirath_matrix_ff_horizontal_concat(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
                                        const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2);

/* Split 'matrix' as 'matrix = [matrix1 | matrix2]. */
void mirath_matrix_ff_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
                                       const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2);

#endif
