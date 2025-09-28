#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include "data_type.h"

/// GF(2^12) modulus
#define MODULUS 0x1009

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
    for (uint32_t j = 0; j < n; j++) {
        // this works only for (q=2, mu=12) and (q=16, mu=3)
        v[j] &= (ff_mu_t)0x0FFF;
    }
}

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    ff_mu_t result = -(a & 1) & b;
    ff_mu_t tmp = b;
    for(int i=1 ; i<12 ; i++) {
        tmp = ((tmp << 1) ^ (-(tmp >> 11) & MODULUS));
        result = result ^ (-(a >> i & 1) & tmp);
    }
    return result;
}

/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    ff_mu_t result = a;
    for(int i=0 ; i<10 ; i++) {
        result = mirath_ff_mu_mult(result, result);
        result = mirath_ff_mu_mult(result, a);
    }
    result = mirath_ff_mu_mult(result, result);
    return result;
}

#endif
