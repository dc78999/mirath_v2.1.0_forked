#ifndef VECTOR_FF_MU_12_H
#define VECTOR_FF_MU_12_H

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
 * \param[out] vector1 out vector1 Vector over ff_mu
 * \param[in] vector2 Vector over ff_mu
 * \param[in] vector3 Vector over ff (in packed representation)
 * \param[in] ncols number of columns
 */
static inline void mirath_vector_ff_mu_add_ff(ff_mu_t *vector1, const ff_mu_t *vector2, const ff_t *vector3, const uint32_t ncols) {

    uint32_t i = ncols;
    const uint32_t limit = ncols % 8;

    while (i >= 16u) {
        const uint8_t t11 = *(vector3 + 0);
        const uint8_t t12 = *(vector3 + 1);

        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);

        const __m128i t1 = _mm_set_epi64x(t22, t21);
        const __m256i m2 = _mm256_cvtepu8_epi16(t1);
        const __m256i m1 = _mm256_loadu_si256((const __m256i *)vector2);

        _mm256_storeu_si256((__m256i *)vector1, m1 ^ m2);
        i       -= 16u;
        vector3 += 2u;
        vector2 += 16u;
        vector1 += 16u;
    }

    while (i >= 8u) {
        const uint32_t t11 = *vector3;
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i m2 = _mm_cvtepi8_epi16(t1);
        const __m128i m1 = _mm_loadu_si128((__m128i *)vector2);

        _mm_storeu_si128((__m128i *)vector1, m1 ^ m2);
        i       -= 8u;
        vector3 += 1u;
        vector2 += 8u;
        vector1 += 8u;
    }

    if (limit) {
        uint16_t tmp[8] __attribute__((aligned(16))) = {0};
        for (uint32_t j = 0; j < limit; j++) { tmp[j] = vector2[j];}
        const __m128i m1 = _mm_load_si128((__m128i *)tmp);

        uint64_t mask = (1ul << (i * 8)) - 1ul;
        const uint32_t t11 = vector3[0];
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101&mask);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i m2 = _mm_cvtepi8_epi16(t1);

        _mm_store_si128((__m128i *)tmp, m2 ^ m1);

        for (uint32_t j = 0; j < limit; j++) { vector1[j] = tmp[j];}
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
    const __m256i A256 = _mm256_set1_epi16(scalar);
    const __m128i A128 = _mm_set1_epi16(scalar);

    while (i >= 16u) {
        const uint8_t t11 = *(vector2 + 0);
        const uint8_t t12 = *(vector2 + 1);

        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);

        const __m128i t1 = _mm_set_epi64x(t22, t21);
        const __m256i mi = _mm256_cvtepu8_epi16(t1);
        const __m256i t2 = mirath_ff_mu_mul_ff_u256(A256, mi);

        _mm256_storeu_si256((__m256i *)vector1, t2);
        i       -= 16u;
        vector2 += 2u;
        vector1 += 16u;
    }

    while (i >= 8u) {
        const uint32_t t11 = *vector2;
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i mi = _mm_cvtepi8_epi16(t1);
        const __m128i t2 = mirath_ff_mu_mul_ff_u128(A128, mi);

        _mm_storeu_si128((__m128i *)vector1, t2);
        i       -= 8u;
        vector2 += 1u;
        vector1 += 8u;
    }

    if (i) {
        ff_mu_t tmp[8] __attribute__((aligned(16))) = {0};

        uint64_t mask = (1ul << (i * 8)) - 1ul;
        const uint32_t t11 = vector2[0];
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101&mask);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i mi = _mm_cvtepi8_epi16(t1);
        const __m128i t2 = mirath_ff_mu_mul_ff_u128(A128, mi);

        _mm_store_si128((__m128i *)tmp, t2);

        for (uint32_t j = 0; j < i; j++) { vector1[j] = tmp[j];}
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
    const uint32_t limit = ncols % 8;
    const __m256i A256 = _mm256_set1_epi16(scalar);
    const __m128i A128 = _mm_set1_epi16(scalar);

    while (i >= 16u) {
        const uint8_t t11 = *(vector3 + 0);
        const uint8_t t12 = *(vector3 + 1);
        const __m256i m = _mm256_loadu_si256((const __m256i *)vector2);

        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);

        const __m128i t1 = _mm_set_epi64x(t22, t21);
        const __m256i mi = _mm256_cvtepu8_epi16(t1);
        const __m256i t2 = mirath_ff_mu_mul_ff_u256(A256, mi);

        _mm256_storeu_si256((__m256i *)vector1, t2^m);
        i       -= 16u;
        vector3 += 2u;
        vector2 += 16u;
        vector1 += 16u;
    }

    while (i >= 8u) {
        const uint32_t t11 = *vector3;
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const __m128i m  = _mm_loadu_si128((const __m128i *)vector2);

        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i mi = _mm_cvtepi8_epi16(t1);
        const __m128i t2 = mirath_ff_mu_mul_ff_u128(A128, mi);

        _mm_storeu_si128((__m128i *)vector1, t2^m);
        i       -= 8u;
        vector3 += 1u;
        vector2 += 8u;
        vector1 += 8u;
    }

    if (limit) {
        ff_mu_t tmp[16] __attribute__((aligned(16))) = {0};
        for (uint32_t j = 0; j < limit; j++) { tmp[j] = vector2[j];}
        const __m128i m  = _mm_load_si128((const __m128i *)tmp);

        uint64_t mask = (1ul << (i * 8)) - 1ul;
        const uint32_t t11 = vector3[0];
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101&mask);
        const __m128i t1 = _mm_set_epi64x(0, t21);
        const __m128i mi = _mm_cvtepi8_epi16(t1);
        const __m128i t2 = mirath_ff_mu_mul_ff_u128(A128, mi);

        _mm_store_si128((__m128i *)tmp, t2^m);

        for (uint32_t j = 0; j < limit; j++) { vector1[j] = tmp[j];}
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

    if (i) {
        uint16_t tmp[32] __attribute__((aligned(16)));
        for (uint32_t j = 0; j < i; ++j) { tmp[j] = vector2[j]; }
        const __m256i m = _mm256_load_si256((const __m256i *)tmp);
        const __m256i t = mirath_ff_mu_mul_u256(m, a256);
        _mm256_store_si256((__m256i *)tmp, t);
        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp[j]; }
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

    // avx2 code
    while (i >= 16u) {
        const __m256i t2 = _mm256_loadu_si256((__m256i *)vector3);
        const __m256i t1 = _mm256_loadu_si256((__m256i *)vector2);
        const __m256i t3 = mirath_ff_mu_mul_u256(t2, a256);
        _mm256_storeu_si256((__m256i *)vector1, t3^t1);
        i       -= 16u;
        vector3 += 16u;
        vector2 += 16u;
        vector1 += 16u;
    }

    if (i) {
        uint16_t tmp[32] __attribute__((aligned(32)));
        for (uint32_t j = 0; j < i; ++j) { tmp[j] = vector3[j]; }
        const __m256i m1 = _mm256_load_si256((const __m256i *)tmp);

        for (uint32_t j = 0; j < i; ++j) { tmp[j] = vector2[j]; }
        const __m256i m2 = _mm256_load_si256((const __m256i *)tmp);

        const __m256i t = mirath_ff_mu_mul_u256(m1, a256);
        _mm256_store_si256((__m256i *)tmp, t^m2);
        for (uint32_t j = 0; j < i; ++j) { vector1[j] = tmp[j]; }
    }
}

///
/// \param a
/// \param b
/// \param n
/// \return
/*static inline ff_mu_t mirath_vector_ff_mu_mul_acc(const ff_mu_t *a, const ff_mu_t *b, const uint32_t n) {
    uint32_t i = n;

    __m256i acc = _mm256_setzero_si256();
    // avx2 code
    while (i >= 16u) {
        const __m256i ta = _mm256_loadu_si256((__m256i *)a);
        const __m256i tb = _mm256_loadu_si256((__m256i *)b);
        const __m256i t1 = mirath_ff_mu_mul_u256(ta, tb);
        acc ^= t1;

        i -= 16u;
        a += 16u;
        b += 16u;
    }


    if (i) {
        uint16_t tmpa[16] __attribute__((aligned(32))) = {0},
                 tmpb[16] __attribute__((aligned(32))) = {0};
        for (uint32_t j = 0; j < i; j++) {
            tmpa[j] = a[j];
            tmpb[j] = b[j];
        }

        const __m256i ta = _mm256_load_si256((__m256i *)tmpa);
        const __m256i tb = _mm256_load_si256((__m256i *)tmpb);
        const __m256i t1 = mirath_ff_mu_mul_u256(ta, tb);
        acc ^= t1;
    }

    const ff_mu_t t = mirath_ff_mu_hadd_u256(acc);
    return t;
}*/

/// bit extends the input vector from gf2 to gf2to12
/// \param out[out]: output vector
/// \param in[in]: input vector
/// \param n number of elements in in/out
static inline void mirath_vector_ff_mu_set_to_ff(ff_mu_t *out,
                                                 const ff_t *in,
                                                 const uint32_t n) {
    uint32_t bytes = n/8;
    const uint32_t limit = n % 16;
    while (bytes >= 2u) {
        const uint8_t t11 = *(in + 0);
        const uint8_t t12 = *(in + 1);
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
        const __m128i t1 = _mm_set_epi64x(t22, t21);
        const __m256i mi = _mm256_cvtepu8_epi16(t1);
        _mm256_storeu_si256((__m256i *)out, mi);
        in  += 2u;
        out += 16u;
        bytes -= 2;
    }
    // bytes += (n & 0x7) ? 1 : 0;
    if (limit) {
        uint16_t tmp[16] __attribute__((aligned(32))) = {0};
        uint8_t t11, t12=0;
        t11 = *(in + 0);
        if(limit > 8) { t12 = *(in + 1); }
        const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
        const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
        const __m128i t1 = _mm_set_epi64x(t22, t21);
        const __m256i mi = _mm256_cvtepu8_epi16(t1);
        _mm256_store_si256((__m256i *)tmp, mi);
        for (uint32_t j = 0; j < limit; j++) { out[j] = tmp[j];}
    }
}

#endif
