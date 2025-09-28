#ifndef VECTOR_FF_MU_8_H
#define VECTOR_FF_MU_8_H

#include <stdint.h>

#include "ff_mu.h"
#include "vector_ff_arith.h"

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
    size_t i = ncols;

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

    for(uint32_t j = 0; j<i; j++) {
        vector1[j] = vector2[j] ^ vector3[j];
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff (in packed representation)
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    while (i >= 32u) {
        const uint32_t t11 = *(vector3 + 0);
        const uint32_t t12 = *(vector3 + 1);
        const uint32_t t13 = *(vector3 + 2);
        const uint32_t t14 = *(vector3 + 3);

        const __m256i m = _mm256_loadu_si256((const __m256i *)vector2);
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
        const uint64_t t23 = _pdep_u64(t13, 0x0101010101010101);
        const uint64_t t24 = _pdep_u64(t14, 0x0101010101010101);

        const __m256i t = _mm256_setr_epi64x(t21, t22, t23, t24);
        _mm256_storeu_si256((__m256i *)vector1, t^m);

        vector3 += 4u;
        vector2 += 32u;
        vector1 += 32u;
        i   -= 32u;
    }

    while (i >= 8u) {
        const uint32_t t1 = *(vector3 + 0);
        const uint64_t m = *(uint64_t *)vector2;
        const uint64_t t2 = _pdep_u64(t1, 0x0101010101010101);
        *(uint64_t *)vector1 = t2 ^ m;

        vector3 += 1u;
        vector2 += 8u;
        vector1 += 8u;
        i   -= 8u;
    }

    if (i) {
        uint8_t tmp[8] = {0};
        const uint8_t mask = (1u << i) - 1u;
        const uint32_t t1 = (*vector3) & mask;

        for (uint64_t j = 0; j < i; ++j) { tmp[j] = vector2[j]; }
        uint64_t m = *(uint64_t *)tmp;

        const uint64_t t2 = _pdep_u64(t1, 0x0101010101010101);
        *(uint64_t *)tmp = t2 ^ m;
        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp[j]; }
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols)
 * \brief vector1 = vector2 * scalar
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector2 Vector over ff (in packed representation)
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_mult_multiple_ff(ff_mu_t *vector1, const ff_mu_t scalar, const ff_t *vector2, const uint32_t ncols) {
    uint32_t i = ncols;
    const __m256i s = _mm256_set1_epi8(scalar);
    const __m256i zero = _mm256_setzero_si256();

    while (i >= 32u) {
        const uint32_t t11 = *(vector2 + 0);
        const uint32_t t12 = *(vector2 + 1);
        const uint32_t t13 = *(vector2 + 2);
        const uint32_t t14 = *(vector2 + 3);

        const uint64_t t21 = _pdep_u64(t11, 0x8080808080808080);
        const uint64_t t22 = _pdep_u64(t12, 0x8080808080808080);
        const uint64_t t23 = _pdep_u64(t13, 0x8080808080808080);
        const uint64_t t24 = _pdep_u64(t14, 0x8080808080808080);

        const __m256i t1 = _mm256_setr_epi64x(t21, t22, t23, t24);
        const __m256i t2 = _mm256_blendv_epi8(zero, s, t1);

        _mm256_storeu_si256((__m256i *)vector1, t2);

        vector2 += 4u;
        vector1 += 32u;
        i       -= 32u;
    }

    if (i) {
        uint64_t tmp[4] __attribute__((aligned(32)));
        uint8_t *tmp2 = (uint8_t *)tmp;
        for (uint32_t j = 0; j < (i+7)/8; ++j) {
            tmp[j] = _pdep_u64(vector2[j], 0x8080808080808080);
        }

        const __m256i t2 = _mm256_blendv_epi8(zero, s, _mm256_load_si256((const __m256i *)(tmp)));
        _mm256_storeu_si256((__m256i *)tmp, t2);
        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp2[j]; }
    }
}

/**
 * \fn static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols)
 * \brief vector1 = vector2 + scalar * vector3
 *
 * \param[out] vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] scalar Scalar over ff_mu
 * \param[in] vector3 Vector over ff (in packed representation)
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_multiple_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_mu_t scalar, const ff_t *vector3, const uint32_t ncols) {
    uint32_t i = ncols;
    const __m256i s = _mm256_set1_epi8(scalar);
    const __m256i zero = _mm256_setzero_si256();

    while (i >= 32u) {
        const uint32_t t11 = *(vector3 + 0);
        const uint32_t t12 = *(vector3 + 1);
        const uint32_t t13 = *(vector3 + 2);
        const uint32_t t14 = *(vector3 + 3);

        const __m256i m = _mm256_loadu_si256((const __m256i *)vector2);
        const uint64_t t21 = _pdep_u64(t11, 0x8080808080808080);
        const uint64_t t22 = _pdep_u64(t12, 0x8080808080808080);
        const uint64_t t23 = _pdep_u64(t13, 0x8080808080808080);
        const uint64_t t24 = _pdep_u64(t14, 0x8080808080808080);

        const __m256i t1 = _mm256_setr_epi64x(t21, t22, t23, t24);
        const __m256i t2 = _mm256_blendv_epi8(zero, s, t1);

        _mm256_storeu_si256((__m256i *)vector1, t2 ^ m);

        vector3 += 4u;
        vector2 += 32u;
        vector1 += 32u;
        i       -= 32u;
    }

    if (i) {
        uint64_t tmp[4] __attribute__((aligned(32)));
        uint8_t *tmp2 = (uint8_t *)tmp;
        for (uint32_t j = 0; j < (i+7)/8; ++j) {
            tmp[j] = _pdep_u64(vector3[j], 0x8080808080808080);
        }

        const __m256i t2 = _mm256_blendv_epi8(zero, s, _mm256_load_si256((const __m256i *)(tmp)));

        for (uint32_t j = 0; j < i; ++j) { tmp2[j] = vector2[j]; }
        const __m256i m = _mm256_load_si256((const __m256i *)tmp);

        _mm256_storeu_si256((__m256i *)tmp, t2^m);
        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp2[j]; }
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
        uint8_t tmp[16] __attribute__((aligned(32)));
        for (uint32_t k = 0; k < i; k++) { tmp[k] = vector3[k]; }
        const __m128i a = _mm_load_si128((const __m128i *) tmp);
        for (uint32_t k = 0; k < i; k++) { tmp[k] = vector2[k]; }
        const __m128i m = _mm_load_si128((const __m128i *) tmp);
        const __m128i c = mirath_ff_mu_linear_transform_8x8_128b(ml128, mh128, a, mask128);

        _mm_store_si128((__m128i *) tmp, c^m);
        for (uint32_t k = 0; k < i; k++) { vector1[k] = tmp[k]; }
    }
}

/// \param n = number of elements in `in`, not bytes.
static inline void mirath_vector_ff_mu_set_to_ff_u256(ff_mu_t *out,
                                                      const ff_t *in,
                                                      const uint32_t n) {
    uint32_t i = n;
    while (i >= 32u) {
        const uint32_t t11 = *(in + 0);
        const uint32_t t12 = *(in + 1);
        const uint32_t t13 = *(in + 2);
        const uint32_t t14 = *(in + 3);

        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
        const uint64_t t23 = _pdep_u64(t13, 0x0101010101010101);
        const uint64_t t24 = _pdep_u64(t14, 0x0101010101010101);

        const __m256i t = _mm256_setr_epi64x(t21, t22, t23, t24);
        _mm256_storeu_si256((__m256i *)out, t);

        in  += 4u;
        out += 32u;
        i   -= 32u;
    }

    while (i >= 8u) {
        const uint32_t t1 = *(in + 0);
        const uint64_t t2 = _pdep_u64(t1, 0x0101010101010101);
        *(uint64_t *)out = t2;

        in  += 1u;
        out += 8u;
        i   -= 8u;
    }

    if (i) {
        uint8_t tmp[8];
        const uint8_t mask = (1u << i) - 1u;
        const uint32_t t1 = (*in) & mask;
        const uint64_t t2 = _pdep_u64(t1, 0x0101010101010101);
        *(uint64_t *)tmp = t2;
        for (uint32_t j = 0; j < i; j++) {
            out[j] = tmp[j];
        }
    }
}


#endif
