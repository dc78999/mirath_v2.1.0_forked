#ifndef ARITH_FF_H
#define ARITH_FF_H

#include <stdint.h>
#include <immintrin.h>

#include "data_type.h"

/// NOTE: assumes a mod 2
/// \param a input element
/// \return a**{-1}
static inline ff_t mirath_ff_inv(const ff_t a) {
   return a;
}

/// NOTE: assumes that a and b are % 2.
/// \param a input element
/// \param b input element
/// \return a*b % 2
static inline ff_t mirath_ff_product(const ff_t a, const ff_t b) {
   return a & b;
}

/// horizontal sum over 256 elements
/// \param in input elements
/// \return popcnt(in) & 1
static inline ff_t gf2_hadd_u256(const __m256i in) {
    __m256i ret = in;
    ret = ret ^ _mm256_permute2x128_si256(ret, ret, 129); // 0b10000001
    ret = ret ^ _mm256_srli_si256(ret, 8);
    const uint64_t r1 = _mm256_extract_epi64(ret, 0);
    const uint64_t r2 = __builtin_popcountll(r1);
    return r2&1;
}
#endif
