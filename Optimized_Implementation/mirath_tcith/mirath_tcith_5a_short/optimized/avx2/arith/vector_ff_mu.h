#ifndef VECTOR_FF_MU_SHORT_H
#define VECTOR_FF_MU_SHORT_H

#include <stdint.h>
#include "matrix_ff_arith.h"
#include "vector_ff_arith.h"
#include "ff_mu.h"


/**
 * \fn static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    // avx2 code
    while (i >= 16u) {
        _mm256_storeu_si256((__m256i *)vector1,
                            _mm256_loadu_si256((__m256i *)vector2) ^
                            _mm256_loadu_si256((__m256i *)vector3));
        i       -= 16u;
        vector1 += 16u;
        vector2 += 16u;
        vector3 += 16u;
    }

    // sse code
    while(i >= 8u) {
        _mm_storeu_si128((__m128i *)vector1,
                         _mm_loadu_si128((__m128i *)vector2) ^
                         _mm_loadu_si128((__m128i *)vector3));
        i       -= 8u;
        vector1 += 8u;
        vector2 += 8u;
        vector3 += 8u;
    }

    for (; i > 0; --i) {
        *vector1++ = *vector2++ ^ *vector3++;
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    // avx2 code
    while (i >= 16u) {
        const __m256i m2 = mirath_ff_mu_extend_gf16_x32(vector3);
        const __m256i m1 = _mm256_loadu_si256((const __m256i *)vector2);

        _mm256_storeu_si256((__m256i *)vector1, m1 ^ m2);
        i       -= 16u;
        vector3 += 8u;
        vector2 += 16u;
        vector1 += 16u;
    }

    while (i >= 8u) {
        const uint32_t t11 = *((uint32_t *)(vector3 + 0));
        const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i m2 = _mm_cvtepi8_epi16(t1);
        const __m128i m1 = _mm_loadu_si128((__m128i *)vector2);

        _mm_storeu_si128((__m128i *)vector1, m1 ^ m2);

        i       -= 8u;
        vector3 += 4u;
        vector2 += 8u;
        vector1 += 8u;
    }

    if (i < 8) { // to avoid this: vector_ff_mu.h: warning: iteration 8 invokes undefined behavior [-Waggressive-loop-optimizations]
        uint16_t tmp[8] __attribute__((aligned(32))) = {0};
        uint8_t *tmp2 = (uint8_t *)tmp;
        uint32_t *tmp3 = (uint32_t *)tmp;

        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp2[j] = vector3[j]; }

        const uint32_t t11 = *tmp3;
        const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i m2 = _mm_cvtepi8_epi16(t1);
        _mm_store_si128((__m128i *)tmp, m2);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = vector2[j] ^ tmp[j]; }
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols)
 * \brief vector1 = vector2 * scalar
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector2 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols) {
    const __m256i s = _mm256_set1_epi16(scalar);
    uint32_t i = ncols;

    // avx2 code
    while (i >= 16u) {
        const __m256i t1 = mirath_ff_mu_extend_gf16_x32(vector2);
        const __m256i t2 = mirath_ff_mu_mul_u256(t1, s);

        _mm256_storeu_si256((__m256i *)vector1, t2);
        i       -= 16u;
        vector2 += 8u;
        vector1 += 16u;
    }

    if (i) {
        uint8_t tmp[32] __attribute__((aligned(32))) = {0};
        uint16_t *tmp2 = (uint16_t *)tmp;
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp[j] = vector2[j]; }

        const __m256i t1 = mirath_ff_mu_extend_gf16_x32(tmp);
        const __m256i t2 = mirath_ff_mu_mul_u256(t1, s);
        _mm256_store_si256((__m256i *)tmp, t2);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp2[j]; }
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + scalar * vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector3 Vector over ff
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols) {
    const __m256i s = _mm256_set1_epi16(scalar);
    uint32_t i = ncols;

    // avx2 code
    while (i >= 16u) {
        const __m256i t1 = mirath_ff_mu_extend_gf16_x32(vector3);
        const __m256i m1 = _mm256_loadu_si256((const __m256i *)vector2);
        const __m256i t2 = mirath_ff_mu_mul_ff_u256(s, t1);

        _mm256_storeu_si256((__m256i *)vector1, t2^m1);
        i       -= 16u;
        vector3 += 8u;
        vector2 += 16u;
        vector1 += 16u;
    }

    if (i) {
        uint8_t tmp[32] __attribute__((aligned(32))) = {0};
        uint8_t tmp3[32] __attribute__((aligned(32))) = {0};
        uint16_t *tmp2 = (uint16_t *)tmp;
        uint16_t *tmp32 = (uint16_t *)tmp3;
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp[j] = vector3[j]; }

        const __m256i t1 = mirath_ff_mu_extend_gf16_x32(tmp);

        for (uint32_t j = 0; j < i; ++j) { tmp32[j] = vector2[j]; }
        const __m256i m1 = _mm256_load_si256((const __m256i *)tmp32);

        const __m256i t2 = mirath_ff_mu_mul_ff_u256(s, t1);
        _mm256_store_si256((__m256i *)tmp, t2^m1);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp2[j]; }
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const uint32_t ncols)
 * \brief vector1 = scalar * vector2
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_mult_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const uint32_t ncols) {
    uint32_t i = ncols;
    const __m256i a256 = _mm256_set1_epi16(scalar);

    // avx2 code
    while (i >= 16u) {
        const __m256i t2 = _mm256_loadu_si256((__m256i *)vector2);
        const __m256i t3 = mirath_ff_mu_mul_u256(t2, a256);
        _mm256_storeu_si256((__m256i *)vector1, t3);
        i       -= 16u;
        vector2 += 16u;
        vector1 += 16u;
    }

    for (; i > 0; --i) {
        *vector1++ = mirath_ff_mu_mult(scalar, *vector2++);
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t scalar, const ff_mu_t *vector2, const ff_mu_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + scalar * vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector3 Vector over ff_mu
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_multiple(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_mu_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    const __m256i a256 = _mm256_set1_epi16(scalar);
    const uint32_t tail = ncols % 16;

    // avx2 code
    while (i >= 16u) {
        const __m256i t1 = _mm256_loadu_si256((__m256i *)vector2);
        const __m256i t2 = _mm256_loadu_si256((__m256i *)vector3);
        const __m256i t3 = mirath_ff_mu_mul_u256(t2, a256);
        const __m256i t = t1 ^ t3;

        _mm256_storeu_si256((__m256i *)vector1, t);
        i       -= 16u;
        vector1 += 16u;
        vector2 += 16u;
        vector3 += 16u;
    }

    if (tail) {
        uint16_t tmp1[16] __attribute__((aligned(32)));
        uint16_t tmp2[16] __attribute__((aligned(32)));
        for (uint32_t j = 0; j < tail; ++j) { tmp1[j] = vector3[j]; }
        const __m256i t3 = _mm256_load_si256((__m256i *)tmp1);

        for (uint32_t j = 0; j < tail; ++j) { tmp2[j] = vector2[j]; }
        const __m256i t2 = _mm256_load_si256((__m256i *)tmp2);

        const __m256i t1 = mirath_ff_mu_mul_u256(t3, a256);
        const __m256i tt = t2 ^ t1;

        _mm256_store_si256((__m256i *)tmp1, tt);
        for (uint32_t j = 0; j < tail; ++j) { vector1[j] = tmp1[j]; }
    }
}

#endif
