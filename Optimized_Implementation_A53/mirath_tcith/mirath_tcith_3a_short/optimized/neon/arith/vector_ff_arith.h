#ifndef VECTOR_ARITH_FF_H
#define VECTOR_ARITH_FF_H

#include <stdint.h>
#include <stdlib.h>
#include <arm_neon.h>
#include "ff.h"

/// \param arg1[out] arg2 + arg3
/// \param arg2[in] first input vector
/// \param arg3[in] second input vector
/// \param d[in] size of each vector
static inline void mirath_vec_ff_add_arith(ff_t *arg1, const ff_t *arg2, const ff_t *arg3, const uint32_t d) {
    uint32_t i = d;
    // neon code
    while(i >= 16u) {
        vst1q_u8((uint8_t *)arg1,
                 veorq_u8(vld1q_u8((const uint8_t *)arg2),
                          vld1q_u8((const uint8_t *)arg3)));
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
