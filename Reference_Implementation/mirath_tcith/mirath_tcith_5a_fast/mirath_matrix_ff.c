#include <stdint.h>
#include <string.h>
#include "mirath_matrix_ff.h"
#include "mirath_parsing.h"

void mirath_matrix_init_zero(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    memset(matrix, 0, mirath_matrix_ff_bytes_size(n_rows, n_cols));
}

void mirath_matrix_expand_seed_public_matrix(ff_t H[MIRATH_VAR_FF_H_BYTES], const seed_t seed_pk) {
    shake_prng_t seedexpander_shake;
    seedexpander_shake_init(&seedexpander_shake, seed_pk, MIRATH_SECURITY_BYTES, NULL, 0);
    seedexpander_shake_get_bytes(&seedexpander_shake, H, MIRATH_VAR_FF_H_BYTES);
    mirath_matrix_set_to_ff(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K);
}

void mirath_matrix_expand_seed_secret_matrix(ff_t S[MIRATH_VAR_FF_S_BYTES], ff_t C[MIRATH_VAR_FF_C_BYTES], const seed_t seed_sk){
    shake_prng_t seedexpander_shake;
    seedexpander_shake_init(&seedexpander_shake, seed_sk, MIRATH_SECURITY_BYTES, NULL, 0);
    seedexpander_shake_get_bytes(&seedexpander_shake, S, MIRATH_VAR_FF_S_BYTES);
    mirath_matrix_set_to_ff(S, MIRATH_PARAM_M, MIRATH_PARAM_R);
    seedexpander_shake_get_bytes(&seedexpander_shake, C, MIRATH_VAR_FF_C_BYTES);
    mirath_matrix_set_to_ff(C, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
}

void mirath_matrix_compute_y(ff_t y[MIRATH_VAR_FF_Y_BYTES],
                                    const ff_t S[MIRATH_VAR_FF_S_BYTES],
                                    const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                    const ff_t H[MIRATH_VAR_FF_H_BYTES]) {

    ff_t e_A[MIRATH_VAR_FF_Y_BYTES] = {0};
    ff_t e_B[mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1)] = {0};

    ff_t T[MIRATH_VAR_FF_T_BYTES] = {0};
    ff_t E[MIRATH_VAR_FF_E_BYTES] = {0};

    mirath_matrix_ff_product(T, S, C, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);
    mirath_matrix_ff_horizontal_concat(E, S, T, MIRATH_PARAM_M, MIRATH_PARAM_R, MIRATH_PARAM_N - MIRATH_PARAM_R);

    const uint32_t bytes_e_B = mirath_matrix_ff_bytes_size(MIRATH_PARAM_K, 1);

    memcpy(e_A, E, MIRATH_VAR_FF_Y_BYTES);
#if (OFF_E_A > 0)
    const uint8_t mask = (1 << (8 - OFF_E_A)) - 1;
    e_A[MIRATH_VAR_FF_Y_BYTES - 1] &= mask;

    for (uint32_t i = 0; i < bytes_e_B - 1 ; i++) {
        e_B[i] = ((E[MIRATH_VAR_FF_Y_BYTES - 1 + i]) >> (8 - OFF_E_A));
        e_B[i] ^= ((E[MIRATH_VAR_FF_Y_BYTES + i]) << (OFF_E_A));
    }
#if ((OFF_E_A + OFF_E_B) >= 8)
    e_B[bytes_e_B - 1] = ((E[MIRATH_VAR_FF_E_BYTES - 1]) >> (8 - OFF_E_A));
#else
    e_B[bytes_e_B - 1] = (E[MIRATH_VAR_FF_E_BYTES - 2] >> (8 - OFF_E_A));
    e_B[bytes_e_B - 1] ^= (E[MIRATH_VAR_FF_E_BYTES - 1] << OFF_E_A);
#endif
#else
    memcpy(e_B, E + MIRATH_VAR_FF_Y_BYTES, bytes_e_B);
#endif

    memset(y, 0, MIRATH_VAR_FF_Y_BYTES);
    mirath_matrix_ff_product(y, H, e_B, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K, 1);

    mirath_vec_ff_add_arith(y, y, e_A, MIRATH_VAR_FF_Y_BYTES);
}

void mirath_tciht_compute_public_key(uint8_t *pk, const uint8_t *sk,
                                     const ff_t S[MIRATH_VAR_FF_S_BYTES],
                                     const ff_t C[MIRATH_VAR_FF_C_BYTES],
                                     const ff_t H[MIRATH_VAR_FF_H_BYTES]) {

    ff_t y[MIRATH_VAR_FF_Y_BYTES];

    mirath_matrix_compute_y(y, S, C, H);

    unparse_public_key(pk, sk + MIRATH_SECURITY_BYTES, y);
}

void mirath_matrix_decompress_secret_key(ff_t S[MIRATH_VAR_FF_S_BYTES],
                      ff_t C[MIRATH_VAR_FF_C_BYTES],
                      ff_t H[MIRATH_VAR_FF_H_BYTES],
                      uint8_t *pk,
                      const uint8_t *sk) {
    seed_t seed_sk;
    seed_t seed_pk;
    ff_t y[MIRATH_VAR_FF_Y_BYTES];

    parse_secret_key(seed_sk, seed_pk, sk);
    mirath_matrix_expand_seed_public_matrix(H, seed_pk);
    mirath_matrix_expand_seed_secret_matrix(S, C,seed_sk);
    mirath_matrix_compute_y(y, S, C, H);
    unparse_public_key(pk, seed_pk, y);
}

void mirath_matrix_decompress_pk(ff_t H[MIRATH_VAR_FF_H_BYTES], ff_t y[MIRATH_VAR_FF_Y_BYTES], const uint8_t *pk) {
    seed_t seed_pk;
    parse_public_key(seed_pk, y, pk);

    shake_prng_t seedexpander_shake;
    seedexpander_shake_init(&seedexpander_shake, (uint8_t*)seed_pk, MIRATH_SECURITY_BYTES, NULL, 0);
    seedexpander_shake_get_bytes(&seedexpander_shake, H, MIRATH_VAR_FF_H_BYTES);
    mirath_matrix_set_to_ff(H, MIRATH_PARAM_M * MIRATH_PARAM_N - MIRATH_PARAM_K, MIRATH_PARAM_K);
}

void mirath_matrix_ff_copy(ff_t *matrix1, const ff_t *matrix2, const uint32_t n_rows, const uint32_t n_cols) {
    const uint32_t n_bytes = mirath_matrix_ff_bytes_size(n_rows, n_cols);
    memcpy(matrix1, matrix2, n_bytes);
}

void mirath_matrix_ff_neg(ff_t *matrix, const uint32_t n_rows, const uint32_t n_cols) {
    /* Nothing to do in characteristic 2. */

    /* Suppress 'unused parameter' warnings. */
    (void)(matrix); (void)(n_rows); (void)(n_cols);
}

void mirath_matrix_ff_add(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols){
    mirath_matrix_ff_add_arith(matrix1, matrix2, matrix3, n_rows, n_cols);
}

void mirath_matrix_ff_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    mirath_matrix_ff_add_multiple_arith(matrix1, scalar, matrix2, n_rows, n_cols);
}

void mirath_matrix_ff_sub(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3, uint32_t n_rows, uint32_t n_cols){
    mirath_matrix_ff_add_arith(matrix1, matrix2, matrix3, n_rows, n_cols);
}

void mirath_matrix_ff_sub_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    mirath_matrix_ff_add_multiple(matrix1, scalar, matrix2, n_rows, n_cols);
}

void mirath_matrix_ff_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    mirath_matrix_ff_product_arith(result, matrix1, matrix2, n_rows1, n_cols1, n_cols2);
}

void mirath_matrix_ff_horizontal_concat(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
                                        const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2) {
    uint8_t *ptr;
    ptr = (uint8_t *)result;
    uint32_t off_ptr = 8;

    uint32_t n_rows_bytes = mirath_matrix_ff_bytes_size(n_rows, 1);
    const uint32_t on_col = 8 - ((8 * n_rows_bytes) - (4 * n_rows));

    uint8_t *col;

    col = (uint8_t *)matrix1;
    for (uint32_t j = 0; j < n_cols1; j++) {
        *ptr |= (*col << (8 - off_ptr));

        for (uint32_t i = 0; i < n_rows_bytes-1; i++) {
            ptr += 1;
            *ptr = (*col >> off_ptr);
            col += 1;
            *ptr |= (*col << (8 - off_ptr));
        }

        if (off_ptr <= on_col) {
            ptr += 1;
            *ptr = (*col >> off_ptr);
        }
        col += 1;
        off_ptr = 8 - ((on_col - off_ptr) % 8);
    }

    col = (uint8_t *)matrix2;
    for (uint32_t j = 0; j < n_cols2; j++) {
        *ptr |= (*col << (8 - off_ptr));

        for (uint32_t i = 0; i < n_rows_bytes-1; i++) {
            ptr += 1;
            *ptr = (*col >> off_ptr);
            col += 1;
            *ptr |= (*col << (8 - off_ptr));
        }

        if (off_ptr <= on_col) {
            ptr += 1;
            if (off_ptr < on_col) {
                *ptr = (*col >> off_ptr);
            }
        }
        col += 1;
        off_ptr = 8 - ((on_col - off_ptr) % 8);
    }
}

void mirath_matrix_ff_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
    const uint32_t n_rows, const uint32_t n_cols1, const uint32_t n_cols2) {
    const uint32_t n_bytes1 = mirath_matrix_ff_bytes_size(n_rows, n_cols1);
    const uint32_t n_bytes2 = mirath_matrix_ff_bytes_size(n_rows, n_cols2);

    if (matrix1 != NULL) {
        memcpy(matrix1, matrix, n_bytes1);
        memcpy(matrix2, matrix + n_bytes1, n_bytes2);
    }
    else {
        memcpy(matrix2, matrix + n_bytes1, n_bytes2);
    }
}
