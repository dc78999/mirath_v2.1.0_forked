#ifndef VECTOR_ARITH_FF_H
#define VECTOR_ARITH_FF_H

#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>

#include "ff.h"

/// NOTE: no boundary checks are performed
/// \param vector[in]: intput vector 
/// \param i[in]: index to return
/// \return: vector[i]
static inline ff_t mirath_vec_ff_get_entry(const ff_t *vector, const uint32_t i) {
	return (vector[i / 8] >> (i % 8)) & 0x01;
}

/// NOTE: no boundary checks are performed.
///     vector[i] = scalar
/// \param vector[in/out]: input vector to alter
/// \param i[in]: position in the vector to alter 
/// \param scalar[in]: value to write
static inline void mirath_vec_ff_set_entry(ff_t *vector, const uint32_t i, const ff_t scalar) {
    const uint8_t mask = 0xff ^ (1 << (i % 8));
    vector[i/8] = (vector[i/8] & mask) ^ (scalar << (i % 8));
}

/// out = in1 + in2
/// \param out[out]: out vector
/// \param in1[in]: input vector 1 
/// \param in2[in]: input vector 2
/// \param d[in]: number of bytes in each vector
static inline void mirath_vec_ff_add_arith(ff_t *arg1, const ff_t *arg2, const ff_t *arg3, const uint32_t d) {
    uint32_t nbytes = d;

    // avx2 code
    while (nbytes >= 32u) {
        _mm256_storeu_si256((__m256i *)arg1,
                            _mm256_loadu_si256((__m256i *)arg2) ^
                            _mm256_loadu_si256((__m256i *)arg3));
        nbytes -= 32u;
        arg1 += 32u;
        arg2 += 32u;
        arg3 += 32u;
    }

    // sse code
    while(nbytes >= 16u) {
        _mm_storeu_si128((__m128i *)arg1,
                         _mm_loadu_si128((__m128i *)arg2) ^
                         _mm_loadu_si128((__m128i *)arg3));
        nbytes -= 16u;
        arg1 += 16u;
        arg2 += 16u;
        arg3 += 16u;
    }

    for (; nbytes > 0; --nbytes) {
        *arg1++ = *arg2++ ^ *arg3++;
    }
}

/// arg1 = arg2 +:while scalar*args2
/// \param arg1
/// \param arg2
/// \param scalar
/// \param arg3
/// \param d[in]: number of bytes in each vector
static inline void mirath_vec_ff_add_multiple_arith(ff_t *arg1, const ff_t *arg2, const ff_t scalar, const ff_t *arg3, const uint32_t d) {
    uint32_t nbytes = d -8;
    const uint8_t m = -scalar;
    const __m256i s256 = _mm256_set1_epi8(m);
    const __m128i s128 = _mm_set1_epi8(m);

    // avx2 code
    while (nbytes >= 32u) {
        const __m256i t1 = _mm256_loadu_si256((const __m256i *)arg2);
        const __m256i t2 = _mm256_loadu_si256((const __m256i *)arg3);
        const __m256i t3 = t1 ^ (t2&s256);
        _mm256_storeu_si256((__m256i *)arg1, t3);

        nbytes -= 32u;
        arg1 += 32u;
        arg2 += 32u;
        arg3 += 32u;
    }

    // sse code
    while(nbytes >= 16u) {
        const __m128i t1 = _mm_loadu_si128((const __m128i *)arg2);
        const __m128i t2 = _mm_loadu_si128((const __m128i *)arg3);
        const __m128i t3 = t1 ^ (t2&s128);
        _mm_storeu_si128((__m128i *)arg1, t3);

        nbytes -= 16u;
        arg1 += 16u;
        arg2 += 16u;
        arg3 += 16u;
    }

    for (; nbytes > 0; --nbytes) {
        *arg1++ = *arg2++ ^ (*arg3++ & m);

    }
}

/// inteprets the vector as a polynomial and evalutates it on `value`
/// \param arg1[in]: the polynomial to evalutate
/// \param value[in]: input value
/// \param d[in]: number of bytes in each vector
/// \return arg1(value)
static inline ff_t mirath_vec_ff_eval_arith(const ff_t *arg1, const ff_t value, const uint32_t d) {
    uint32_t nbytes = d;
    const uint8_t m = -value;
    const __m256i s256 = _mm256_set1_epi8(m);

    __m256i acc = _mm256_setzero_si256();
    // avx2 code
    while (nbytes >= 32u) {
        const __m256i t1 = _mm256_loadu_si256((const __m256i *)arg1);
        const __m256i t2 = t1&s256;
        acc ^= t2;

        nbytes -= 32u;
        arg1 += 32u;
    }

    
    if (nbytes) {
        uint8_t tmp[32] __attribute__((aligned(32))) = {0};
        for (uint32_t i = 0; i < nbytes; i++) { tmp[i] = arg1[i]; }

        __m256i t1 = _mm256_load_si256((const __m256i *)tmp);
        const __m256i t2 = t1&s256;
        acc ^= t2;
    }
    
    return gf2_hadd_u256(acc);
}
#endif
