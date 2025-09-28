#ifndef VECTOR_ARITH_FF_H
#define VECTOR_ARITH_FF_H

#include <stdint.h>
#include <stdlib.h>
#include "ff.h"

static inline ff_t mirath_vec_ff_get_entry(const ff_t *vector, const uint32_t i) {
	return (vector[i / 8] >> (i % 8)) & 0x01;
}

static inline void mirath_vec_ff_set_entry(ff_t *vector, const uint32_t i, const ff_t scalar) {
    const uint8_t mask = 0xff ^ (1 << (i % 8));
    vector[i/8] = (vector[i/8] & mask) ^ (scalar << (i % 8));
}

static inline void mirath_vec_ff_add_arith(ff_t *arg1, const ff_t *arg2, const ff_t *arg3, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = arg2[i] ^ arg3[i];
	}
}

static inline void mirath_vec_ff_add_multiple_arith(ff_t *arg1, const ff_t *arg2, const ff_t scalar, const ff_t *arg3, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = arg2[i] ^ (arg3[i] & (-scalar));
	}
}

static inline void mirath_vec_ff_scalar_mult_arith(ff_t *arg1, const ff_t *arg2, const ff_t value, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] = arg2[i] & (-value);
	}
}

static inline void mirath_vec_ff_scalar_add_arith(ff_t *arg1, const ff_t *arg2, const ff_t value, const uint32_t d) {
	for (uint32_t i = 0; i < d; i++) {
		arg1[i] ^= arg2[i] & (-value);
	}
}

static inline ff_t * mirath_vec_ff_mult_arith(ff_t *arg1, const ff_t *arg2, const uint32_t d) {
	ff_t *ret = (ff_t *)calloc(d, 1);
	for (uint32_t i = 0; i < d; i++) {
		ff_t coeff = 0;
		for (uint32_t j = 0; j < d; j++) {
			ret[i] ^= mirath_ff_product(mirath_vec_ff_get_entry(arg1, i), mirath_vec_ff_get_entry(arg2, j));
		}
		mirath_vec_ff_set_entry(ret, i, coeff);
	}

	return ret;
}

static inline ff_t mirath_vec_ff_eval_arith(const ff_t *arg1, const ff_t value, const uint32_t d) {
	ff_t ret = 0;
	for (int i = d-1; i>=0; i--) {
		ret ^= mirath_ff_product(value, mirath_vec_ff_get_entry(arg1, i));
	}

	return ret;
}
#endif
