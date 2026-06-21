#ifndef VECTOR_FF_MU_FAST_H
#define VECTOR_FF_MU_FAST_H

#include <stdint.h>
#include <arm_neon.h>
#include "matrix_ff_arith.h"
#include "vector_ff_arith.h"
#include "ff_mu.h"

static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    while (i >= 16u) {
        vst1q_u8((uint8_t *)vector1,
                 veorq_u8(vld1q_u8((const uint8_t *)vector2),
                          vld1q_u8((const uint8_t *)vector3)));
        i       -= 16u;
        vector1 += 16u;
        vector2 += 16u;
        vector3 += 16u;
    }
    for (; i > 0; --i) {
        *vector1++ = *vector2++ ^ *vector3++;
    }
}

static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    while (i >= 16u) {
        uint8x16_t m = vld1q_u8((const uint8_t *)vector2);
        uint8x16_t t = mirath_ff_mu_extend_gf16_x16(vector3);
        vst1q_u8((uint8_t *)vector1, veorq_u8(t, m));

        vector3 += 8u;
        vector2 += 16u;
        vector1 += 16u;
        i   -= 16u;
    }

    if (i) {
        uint8_t tmp1[16] __attribute__((aligned(16))) = { 0 };
        uint8_t tmp2[16] __attribute__((aligned(16))) = { 0 };
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp1[j] = vector3[j]; }
        for (uint64_t j = 0; j < i; ++j) { tmp2[j] = vector2[j]; }

        uint8x16_t m = vld1q_u8(tmp2);
        uint8x16_t t = mirath_ff_mu_extend_gf16_x16(tmp1);
        vst1q_u8(tmp1, veorq_u8(t, m));
        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp1[j]; }
    }
}

static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols) {
    uint32_t i = ncols;
    uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(scalar);
    uint8x16_t mask = vdupq_n_u8(0xf);

    while (i >= 16u) {
        uint8x16_t a = mirath_ff_mu_extend_gf16_x16(vector2);
        uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8((uint8_t *)vector1, t);

        vector2 += 8u;
        vector1 += 16u;
        i   -= 16u;
    }

    if (i) {
        uint8_t tmp[16] __attribute__((aligned(16))) = { 0 };
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp[j] = vector2[j]; }

        uint8x16_t a = mirath_ff_mu_extend_gf16_x16(tmp);
        uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8(tmp, t);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp[j]; }
    }
}

static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(scalar);
    uint8x16_t mask = vdupq_n_u8(0xf);

    while (i >= 16u) {
        uint8x16_t m = vld1q_u8((const uint8_t *)vector2);
        uint8x16_t a = mirath_ff_mu_extend_gf16_x16(vector3);
        uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8((uint8_t *)vector1, veorq_u8(t, m));

        vector3 += 8u;
        vector2 += 16u;
        vector1 += 16u;
        i   -= 16u;
    }

    if (i) {
        uint8_t tmp1[16] __attribute__((aligned(16))) = {0};
        uint8_t tmp2[16] __attribute__((aligned(16))) = {0};

        for (uint32_t j = 0; j < i; ++j) { tmp1[j] = vector2[j]; }
        uint8x16_t m = vld1q_u8(tmp1);

        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp2[j] = vector3[j]; }
        uint8x16_t a = mirath_ff_mu_extend_gf16_x16(tmp2);
        uint8x16_t t = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8(tmp1, veorq_u8(t, m));

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp1[j]; }
    }
}

static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const uint32_t ncols) {
    size_t i = ncols;
    uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(scalar);
    uint8x16_t mask = vdupq_n_u8(0xf);

    while (i >= 16u) {
        uint8x16_t a = vld1q_u8((const uint8_t *)vector2);
        uint8x16_t tmp = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8((uint8_t *)vector1, tmp);
        vector2 += 16u;
        vector1 += 16u;
        i       -= 16u;
    }

    if (i) {
        uint8_t tmp[16] __attribute__((aligned(16))) = {0};
        for (uint32_t k = 0; k < i; k++) { tmp[k] = vector2[k]; }
        uint8x16_t a = vld1q_u8(tmp);
        uint8x16_t c = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);

        vst1q_u8(tmp, c);
        for (uint32_t k = 0; k < i; k++) { vector1[k] = tmp[k]; }
    }
}

static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t ncols) {
    size_t i = ncols;
    uint8x16x2_t tab = mirath_ff_mu_generate_multab_16_single_element_u128(scalar);
    uint8x16_t mask = vdupq_n_u8(0xf);

    while (i >= 16u) {
        uint8x16_t a = vld1q_u8((const uint8_t *)vector3);
        uint8x16_t tmp = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);
        vst1q_u8((uint8_t *)vector1, veorq_u8(tmp, vld1q_u8((const uint8_t *)vector2)));
        vector3 += 16u;
        vector2 += 16u;
        vector1 += 16u;
        i       -= 16u;
    }

    if (i) {
        uint8_t tmp1[16] __attribute__((aligned(16))) = {0};
        uint8_t tmp2[16] __attribute__((aligned(16))) = {0};
        for (uint32_t k = 0; k < i; k++) { tmp1[k] = vector3[k]; }
        uint8x16_t a = vld1q_u8(tmp1);
        for (uint32_t k = 0; k < i; k++) { tmp2[k] = vector2[k]; }
        uint8x16_t m = vld1q_u8(tmp2);
        uint8x16_t c = mirath_ff_mu_linear_transform_8x8_128b(tab.val[0], tab.val[1], a, mask);

        vst1q_u8(tmp1, veorq_u8(c, m));
        for (uint32_t k = 0; k < i; k++) { vector1[k] = tmp1[k]; }
    }
}

#endif
