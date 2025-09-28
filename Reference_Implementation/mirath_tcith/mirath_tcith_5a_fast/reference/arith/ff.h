#ifndef ARITH_FF_H
#define ARITH_FF_H

#include <stdint.h>
#include "data_type.h"


#define FF_MODULUS 3u

/// i*{-1} = [i]
static const ff_t mirath_ff_inv_table[16] __attribute__((aligned(16))) = {
    0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8
};

static inline ff_t mirath_ff_add(const ff_t a,
                                 const ff_t b) {
    return a ^ b;
}

/// NOTE: assumes a mod 16
/// \param a input element
/// \return a**{-1}
static inline ff_t mirath_ff_inv(const ff_t a) {
    return mirath_ff_inv_table[a];
}

/// NOTE: assumes that a and b are % 16.
/// \param a input element
/// \param b input element
/// \return a*b % 16
static inline ff_t mirath_ff_product(const ff_t a, const ff_t b) {
    uint8_t r;
    r =    (-(b>>3    ) & a);
    r =    (-(b>>2 & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
    r =    (-(b>>1 & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
    return (-(b    & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
}
#endif
