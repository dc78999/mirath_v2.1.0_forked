#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <immintrin.h>
#include "data_type.h"

/// GF(2^12) modulus
#define MODULUS 0x1009

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
    for (uint32_t j = 0; j < n; j++) {
        // this works only for (q=2, mu=12) and (q=16, mu=3)
        v[j] &= (ff_mu_t)0x0FFF;
    }
}

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    ff_mu_t result = -(a & 1) & b;
    ff_mu_t tmp = b;
    for(int i=1 ; i<12 ; i++) {
        tmp = ((tmp << 1) ^ (-(tmp >> 11) & MODULUS));
        result = result ^ (-(a >> i & 1) & tmp);
    }
    return result;
}

/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    ff_mu_t result = a;
    for(int i=0 ; i<10 ; i++) {
        result = mirath_ff_mu_mult(result, result);
        result = mirath_ff_mu_mult(result, a);
    }
    result = mirath_ff_mu_mult(result, result);
    return result;
}

/// \param in[in]: 2bytes which will be extended to a avx2 register
///     = [in0, ..., in_15], where in_i is just a single bit
/// \return [in_0, ..., in_15], extened to 16bit limbs
static inline __m256i mirath_ff_mu_expand_ff_x16_u256(const uint8_t *in) {
    const uint8_t t11 = *(in + 0);
    const uint8_t t12 = *(in + 1);

    const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
    const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);

    const __m128i t1 = _mm_set_epi64x(t22, t21);
    return _mm256_cvtepu8_epi16(t1);
}

/// \param in[in]: bytes which will be extended to a sse register
///     = [in0, ..., in_7], where in_i is just a single bit
/// \return [in_0, ..., in_7], extened to 16bit limbs
static inline __m128i mirath_ff_mu_expand_ff_x8_u256(const uint8_t *in) {
    const uint32_t t11 = *in;
    const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
    const __m128i t1 = _mm_set_epi64x(0, t21);
    return _mm_cvtepi8_epi16(t1);
}


/// horizontal add, but not withing a single limb, but over the 16 -16bit limbs
/// NOTE: add = xor
/// \param in[in]: [in_0, ..., in_15]
/// \return sum(in_0, in_1, ..., in[15])
static inline uint16_t mirath_ff_mu_hadd_u256(const __m256i in) {
    __m256i ret = _mm256_xor_si256(in, _mm256_srli_si256(in, 2));
    ret = _mm256_xor_si256(ret, _mm256_srli_si256(in, 4));
    ret = _mm256_xor_si256(ret, _mm256_srli_si256(ret, 8));
    ret = _mm256_xor_si256(ret, _mm256_permute2x128_si256(ret, ret, 129)); // 0b10000001
    return _mm256_extract_epi16(ret, 0);
}

/// full 16x multiplication
/// \param a[in] = [a_0, ..., a_15]
/// \param b[in] = [b_0, ..., b_15]
/// \return [a_0*b_0, ..., a_15*b_15]
static inline __m256i mirath_ff_mu_mul_u256(const __m256i a,
                                            const __m256i b) {
    const __m256i mod  = _mm256_set1_epi16((short)MODULUS);
    const __m256i one  = _mm256_set1_epi8(-1);
    const __m256i zero = _mm256_set1_epi8(0);
    const __m256i m1 = _mm256_set1_epi16(0x8000);
    const __m256i m2 = _mm256_set1_epi16(0x0080);
    __m256i ma, mr, r;


    // 11
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 4) & m1) ^ (_mm256_srli_epi16(b, 4)));
    r = ma & a;
    // 10
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 5) & m1) ^ (_mm256_srli_epi16(b, 3)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 9
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 6) & m1) ^ (_mm256_srli_epi16(b, 2)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 8
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 7) & m1) ^ (_mm256_srli_epi16(b, 1)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 7
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 8)) ^ (b & m2));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 6
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 9)) ^ (_mm256_slli_epi16(b, 1)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 5
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 10)) ^ (_mm256_slli_epi16(b, 2)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 4
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 11)) ^ (_mm256_slli_epi16(b, 3)));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 3
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 12)) ^ (_mm256_slli_epi16(b, 4) & m2));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 2
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 13)) ^ (_mm256_slli_epi16(b, 5) & m2));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 1
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 14)) ^ (_mm256_slli_epi16(b, 6) & m2));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    // 0
    mr = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(r, 4) & m1) ^ (_mm256_srli_epi16(r, 4)));
    ma = _mm256_blendv_epi8(zero, one, (_mm256_slli_epi16(b, 15)) ^ (_mm256_slli_epi16(b, 7) & m2));
    r =  (ma & a) ^ (mr & mod) ^ _mm256_add_epi16(r, r);
    return r;
}

/// \param a in gf2to12
/// \param b in gf2, not compresses: a single bit in
/// \return a*b
static inline __m256i mirath_ff_mu_mul_ff_u256(const __m256i a,
                                               const __m256i b) {
    const __m256i m1 = _mm256_set1_epi16(-1);
    const __m256i t1 = _mm256_sign_epi16(b, m1);
    return a & t1;
}

/// \param a in gf2to12
/// \param b in gf2, not compresses: a single bit in
/// \return a*b
static inline __m128i mirath_ff_mu_mul_ff_u128(const __m128i a,
                                               const __m128i b) {
    const __m128i m1 = _mm_set1_epi16(-1);
    const __m128i t1 = _mm_sign_epi16(b, m1);
    return a & t1;
}


#endif
