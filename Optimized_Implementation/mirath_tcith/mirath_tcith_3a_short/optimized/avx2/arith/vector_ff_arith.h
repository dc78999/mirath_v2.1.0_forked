#ifndef VECTOR_ARITH_FF_H
#define VECTOR_ARITH_FF_H

#include <stdint.h>
#include <stdlib.h>
#include "ff.h"

/// \param arg1[out] arg2 + arg3
/// \param arg2[in] first input vector
/// \param arg3[in] second input vector
/// \param d[in] size of each vector
static inline void mirath_vec_ff_add_arith(ff_t *arg1, const ff_t *arg2, const ff_t *arg3, const uint32_t d) {
    uint32_t i = d;
    // avx2 code
    while (i >= 32u) {
        _mm256_storeu_si256((__m256i *)arg1,
                            _mm256_loadu_si256((__m256i *)arg2) ^
                            _mm256_loadu_si256((__m256i *)arg3));
        i -= 32u;
        arg1 += 32u;
        arg2 += 32u;
        arg3 += 32u;
    }

    // sse code
    while(i >= 16u) {
        _mm_storeu_si128((__m128i *)arg1,
                         _mm_loadu_si128((__m128i *)arg2) ^
                         _mm_loadu_si128((__m128i *)arg3));
        arg1 += 16u;
        arg2 += 16u;
        arg3 += 16u;
        i -= 16;
    }

    for(uint32_t j = 0; j<i; j++) {
        arg1[j] = arg2[j] ^ arg3[j];
    }
}

#endif
