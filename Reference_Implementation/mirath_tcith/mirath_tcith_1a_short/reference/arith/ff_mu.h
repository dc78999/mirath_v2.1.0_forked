#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <stdint.h>

#include "ff.h"
#include "data_type.h"

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
    for (uint32_t j = 0; j < n; j++) {
        // this works only for (q=2, mu=12) and (q=16, mu=3)
        v[j] &= (ff_mu_t)0x0FFF;
    }
}

static const uint8_t mirath_map_ff_to_ff_mu[16] __attribute__((aligned(32))) = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a,
                                        const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a,
                                         const ff_mu_t b) {
    ff_mu_t r;

    const ff_t a0 = a&0xF, a1 = (a>>4)&0XF, a2 = (a>>8)&0xF;
    const ff_t b0 = b&0xF, b1 = (b>>4)&0XF, b2 = (b>>8)&0xF;

    const ff_t p0 = mirath_ff_product(a0, b0);
    const ff_t p1 = mirath_ff_product(a1, b1);
    const ff_t p2 = mirath_ff_product(a2, b2);

    const ff_t a01 = mirath_ff_add(a0, a1);
    const ff_t a12 = mirath_ff_add(a1, a2);
    const ff_t a02 = mirath_ff_add(a0, a2);
    const ff_t b01 = mirath_ff_add(b0, b1);
    const ff_t b12 = mirath_ff_add(b1, b2);
    const ff_t b02 = mirath_ff_add(b0, b2);
    const ff_t p01 = mirath_ff_product(a01, b01);
    const ff_t p12 = mirath_ff_product(a12, b12);
    const ff_t p02 = mirath_ff_product(a02, b02);

    // compute lowest limb
    r = mirath_ff_add(p1, p2);
    r = mirath_ff_add(r, p12);
    r = mirath_ff_add(r, p0);

    r^= p0 << 4;
    r^= p01 << 4;
    r^= p12 << 4;

    r^= p02 << 8;
    r^= p0 << 8;
    r^= p1 << 8;

    return r;
}


#endif
