#ifndef VECTOR_FF_MU_FAST_H
#define VECTOR_FF_MU_FAST_H

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
    while (i >= 32u) {
        _mm256_storeu_si256((__m256i *)vector1,
                            _mm256_loadu_si256((__m256i *)vector2) ^
                            _mm256_loadu_si256((__m256i *)vector3));
        i       -= 32u;
        vector1 += 32u;
        vector2 += 32u;
        vector3 += 32u;
    }

    // sse code
    while(i >= 16u) {
        _mm_storeu_si128((__m128i *)vector1,
                         _mm_loadu_si128((__m128i *)vector2) ^
                         _mm_loadu_si128((__m128i *)vector3));
        i       -= 16;
        vector1 += 16u;
        vector2 += 16u;
        vector3 += 16u;
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
    const __m128i perm128 = _mm_load_si128((const __m128i *)mirath_map_ff_to_ff_mu);

    uint32_t i = ncols;
    while (i >= 32u) {
        const __m256i m = _mm256_loadu_si256((const __m256i *)vector2);
        const __m256i t = mirath_ff_mu_extend_gf16_x32(vector3);
        _mm256_storeu_si256((__m256i *)vector1, t^m);

        vector3 += 16u;
        vector2 += 32u;
        vector1 += 32u;
        i   -= 32u;
    }

    while (i >= 8u) {
        const uint32_t t1 = *(uint32_t *)(vector3);
        const uint64_t t2 = _pdep_u64(t1, 0x0F0F0F0F0F0F0F0F);

        const uint64_t m = *(uint64_t *)vector2;
        const __m128i t3 = _mm_setr_epi64((__m64)t2, (__m64)0ull);
        const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
        const uint64_t t5 = _mm_extract_epi64(t4, 0);
        *(uint64_t *)vector1 = t5 ^ m;

        vector3 += 4u;
        vector2 += 8u;
        vector1 += 8u;
        i   -= 8u;
    }

    if (i) {
        uint8_t tmp1[8] __attribute__((aligned(32))) = { 0 };
        uint8_t tmp2[8] __attribute__((aligned(32))) = { 0 };
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp1[j] = vector3[j]; }
        const uint32_t t1 = *(uint32_t *)tmp1;

        for (uint64_t j = 0; j < i; ++j) { tmp2[j] = vector2[j]; }
        uint64_t m = *(uint64_t *)tmp2;

        const uint64_t t2 = _pdep_u64(t1, 0x0F0F0F0F0F0F0F0F);
        const __m128i t3 = _mm_setr_epi64((__m64)t2, (__m64)0ull);
        const __m128i t4 = _mm_shuffle_epi8(perm128, t3);
        const uint64_t t5 = _mm_extract_epi64(t4, 0);

        *(uint64_t *)tmp1 = t5 ^ m;
        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp1[j]; }
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
    uint32_t i = ncols;

    const __m256i tab = mirath_ff_mu_generate_multab_16_single_element_u256(scalar);
    const __m256i ml = _mm256_permute2x128_si256(tab, tab, 0);
    const __m256i mh = _mm256_permute2x128_si256(tab, tab, 0x11);
    const __m256i mask = _mm256_set1_epi8(0xf);

    while (i >= 32u) {
        const __m256i a = mirath_ff_mu_extend_gf16_x32(vector2);
        const __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_storeu_si256((__m256i *)vector1, t);

        vector2 += 16u;
        vector1 += 32u;
        i   -= 32u;
    }

    if (i) {
        uint8_t tmp[32] __attribute__((aligned(32)));
        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp[j] = vector2[j]; }

        const __m256i a = mirath_ff_mu_extend_gf16_x32(tmp);
        const __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_store_si256((__m256i *)tmp, t);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp[j]; }
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
    uint32_t i = ncols;

    const __m256i tab = mirath_ff_mu_generate_multab_16_single_element_u256(scalar);
    const __m256i ml = _mm256_permute2x128_si256(tab, tab, 0);
    const __m256i mh = _mm256_permute2x128_si256(tab, tab, 0x11);
    const __m256i mask = _mm256_set1_epi8(0xf);

    while (i >= 32u) {
        const __m256i m = _mm256_loadu_si256((const __m256i *)vector2);
        const __m256i a = mirath_ff_mu_extend_gf16_x32(vector3);
        const __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_storeu_si256((__m256i *)vector1, t ^ m);

        vector3 += 16u;
        vector2 += 32u;
        vector1 += 32u;
        i   -= 32u;
    }

    if (i) {
        uint8_t tmp1[32] __attribute__((aligned(32))) = {0};
        uint8_t tmp2[32] __attribute__((aligned(32))) = {0};

        for (uint32_t j = 0; j < i; ++j) { tmp1[j] = vector2[j]; }
        const __m256i m = _mm256_load_si256((const __m256i *)tmp1);

        for (uint32_t j = 0; j < (i+1)/2; ++j) { tmp2[j] = vector3[j]; }
        const __m256i a = mirath_ff_mu_extend_gf16_x32(tmp2);
        const __m256i t = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_store_si256((__m256i *)tmp1, t ^ m);

        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp1[j]; }
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
    size_t i = ncols;
    const __m256i tab = mirath_ff_mu_generate_multab_16_single_element_u256(scalar);
    const __m256i ml = _mm256_permute2x128_si256(tab, tab, 0);
    const __m256i mh = _mm256_permute2x128_si256(tab, tab, 0x11);
    const __m256i mask = _mm256_set1_epi8(0xf);

    const __m128i ml128 = _mm256_extracti128_si256(ml, 0);
    const __m128i mh128 = _mm256_extracti128_si256(mh, 0);
    const __m128i mask128 = _mm_set1_epi8(0xF);

    // avx2 code
    while (i >= 32u) {
        const __m256i a = _mm256_loadu_si256((__m256i *)vector2);
        const __m256i tmp = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_storeu_si256((__m256i *)vector1, tmp);
        vector2 += 32u;
        vector1 += 32u;
        i       -= 32u;
    }

    // sse code
    while(i >= 16u) {
        const __m128i a = _mm_loadu_si128((__m128i *)vector2);
        const __m128i tmp = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, a, mask128);

        _mm_storeu_si128((__m128i *)vector1, tmp);
        vector2 += 16u;
        vector1 += 16u;
        i       -= 16;
    }

    if (i) {
        uint8_t tmp[16] __attribute__((aligned(32)));
        for (uint32_t k = 0; k < i; k++) { tmp[k] = vector2[k]; }
        const __m128i a = _mm_load_si128((const __m128i *) tmp);
        const __m128i c = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, a, mask128);

        _mm_store_si128((__m128i *) tmp, c);
        for (uint32_t k = 0; k < i; k++) { vector1[k] = tmp[k]; }
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
    size_t i = ncols;
    const __m256i tab = mirath_ff_mu_generate_multab_16_single_element_u256(scalar);
    const __m256i ml = _mm256_permute2x128_si256(tab, tab, 0);
    const __m256i mh = _mm256_permute2x128_si256(tab, tab, 0x11);
    const __m256i mask = _mm256_set1_epi8(0xf);

    const __m128i ml128 = _mm256_extracti128_si256(ml, 0);
    const __m128i mh128 = _mm256_extracti128_si256(mh, 0);
    const __m128i mask128 = _mm_set1_epi8(0xF);

    // avx2 code
    while (i >= 32u) {
        const __m256i a = _mm256_loadu_si256((__m256i *)vector3);
        const __m256i tmp = mirath_ff_mu_linear_transform_8x8_256b(ml, mh, a, mask);
        _mm256_storeu_si256((__m256i *)vector1, tmp ^ _mm256_loadu_si256((__m256i *)vector2));
        vector3 += 32u;
        vector2 += 32u;
        vector1 += 32u;
        i       -= 32u;
    }

    // sse code
    while(i >= 16u) {
        const __m128i a = _mm_loadu_si128((__m128i *)vector3);
        const __m128i tmp = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, a, mask128);

        _mm_storeu_si128((__m128i *)vector1, tmp ^ _mm_loadu_si128((const __m128i *)vector2));
        vector3 += 16u;
        vector2 += 16u;
        vector1 += 16u;
        i       -= 16;
    }

    if (i) {
        uint8_t tmp1[16] __attribute__((aligned(32)));
        uint8_t tmp2[16] __attribute__((aligned(32)));
        for (uint32_t k = 0; k < i; k++) { tmp1[k] = vector3[k]; }
        const __m128i a = _mm_load_si128((const __m128i *) tmp1);
        for (uint32_t k = 0; k < i; k++) { tmp2[k] = vector2[k]; }
        const __m128i m = _mm_load_si128((const __m128i *) tmp2);
        const __m128i c = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, a, mask128);

        _mm_store_si128((__m128i *) tmp1, c^m);
        for (uint32_t k = 0; k < i; k++) { vector1[k] = tmp1[k]; }
    }
}

#endif
