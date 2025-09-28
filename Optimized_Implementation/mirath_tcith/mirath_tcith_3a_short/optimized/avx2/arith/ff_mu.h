#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <stdint.h>

#include "ff.h"
#include "data_type.h"

#include <stdio.h>

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
    for (uint32_t j = 0; j < n; j++) {
        // this works only for (q=2, mu=12) and (q=16, mu=3)
        v[j] &= (ff_mu_t)0x0FFF;
    }
}

static const uint8_t mirath_map_ff_to_ff_mu[16] __attribute__((aligned(32))) = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a,
                                       const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a,
                                        const ff_mu_t b) {
    ff_mu_t r;

    const ff_t a0 = a&0xF, a1 = (a>>4)&0XF, a2 = (a>>8)&0xF;
    const ff_t b0 = b&0xF, b1 = (b>>4)&0XF, b2 = (b>>8)&0xF;

    const ff_t p0 = mirath_ff_product(a0, b0);
    const ff_t p1 = mirath_ff_product(a1, b1);
    const ff_t p2 = mirath_ff_product(a2, b2);

    const ff_t a01 = mirath_ff_add(a0, a1);
    const ff_t a12 = mirath_ff_add(a1, a2);
    const ff_t a02 = mirath_ff_add(a0, a2);
    const ff_t b01 = mirath_ff_add(b0, b1);
    const ff_t b12 = mirath_ff_add(b1, b2);
    const ff_t b02 = mirath_ff_add(b0, b2);
    const ff_t p01 = mirath_ff_product(a01, b01);
    const ff_t p12 = mirath_ff_product(a12, b12);
    const ff_t p02 = mirath_ff_product(a02, b02);

    // compute lowest limb
    r = mirath_ff_add(p1, p2);
    r = mirath_ff_add(r, p12);
    r = mirath_ff_add(r, p0);

    r^= p0 << 4;
    r^= p01 << 4;
    r^= p12 << 4;

    r^= p02 << 8;
    r^= p0 << 8;
    r^= p1 << 8;

    return r;
}

/// loads 16 elements/ 8 bytes and extends them to gf16tp3
static inline __m256i mirath_ff_mu_extend_gf16_x32(const uint8_t *in) {
    const uint32_t t11 = *((uint32_t *)(in + 0));
    const uint32_t t12 = *((uint32_t *)(in + 4));
    const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
    const uint64_t t22 = _pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);
    const __m128i t1 = _mm_set_epi64x(t22, t21);
    return _mm256_cvtepu8_epi16(t1);
}

///
static inline __m256i mirath_ff_mu_mul_u256(const __m256i a,
                                            const __m256i b) {
    __m256i r,t;
    const __m256i m0     = _mm256_set1_epi16(0x00F);
    const __m256i m1     = _mm256_set1_epi16(0x0F0);
    const __m256i m12    = _mm256_set1_epi16(0xFF0);
    const __m256i m      = _mm256_set1_epi16(0xFFF);

    // bit  0     4      8     12     16
    // pi = [a0*b0, a1*b1, a2*b2,  0  ]
    // pi = [  p0 ,   p1 ,   p2 ,  0  ]
    const __m256i pi = mirath_ff_mul_full_u256(a, b);


    // bit    0     4      8     12    16
    // a01_12=[a0^a1, a1^a2, a0^a2,   0 ]  = [a01, a12, a02, 0]
    // b01_12=[b0^b1, b1^b2, a0^b2,   0 ]  = [b01, b12, b02, 0]
    const __m256i a01_12 = a ^ _mm256_srli_epi16(a, 4) ^ _mm256_slli_epi16(a, 8);
    const __m256i b01_12 = b ^ _mm256_srli_epi16(b, 4) ^ _mm256_slli_epi16(b, 8);

    // bit    0   4    8   12    16
    // p012 = [p01, p12, p02,   0]
    __m256i p012 = mirath_ff_mul_full_u256(a01_12, b01_12);


    // bit 0   4           8         12    16
    // r = [  0, p0^p12^p01, p1^p0^p02,   0]
    r  = _mm256_slli_epi16(pi, 4);
    r ^= _mm256_slli_epi16(pi, 8);
    r ^= p012 & m12;
    r ^= _mm256_slli_epi16(p012, 4) & m1;


    t = pi ^ (_mm256_srli_epi16(pi, 4));
    t ^= _mm256_srli_epi16(pi, 8);
    t ^= _mm256_srli_epi16(p012, 4);
    t &= m0;
    r ^= t;
    return r&m;
}

/// \param a[in]: [a_0, ..., a_15], a_i in gf16^3
/// \param b[in]: [b_0, ..., b_15], b_i in gf16
/// \return: [a_0*b_0, ..., a_15*_15]
static inline __m256i mirath_ff_mu_mul_ff_u256(const __m256i a,
                                               const __m256i b) {
    const __m256i c = b ^ _mm256_slli_epi16(b, 8);
    return mirath_ff_mul_u256(a, c);
}

/// loads 8 elements/ 4 bytes and extends them to gf256
static inline __m128i mirath_vector_extend_gf16_x8(const ff_t *in) {
    const uint32_t t11 = *((uint32_t *)(in + 0));
    const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
    const __m128i t1 = _mm_set_epi64x(0, t21);
    return _mm_cvtepi8_epi16(t1);
}

/// loads 16 elements/ 8 bytes and extends them to gf256
static inline __m256i mirath_vector_extend_gf16_x16(const ff_t *in) {
    const uint32_t t11 = *((uint32_t *)(in + 0));
    const uint32_t t12 = *((uint32_t *)(in + 4));
    const uint64_t t21 = _pdep_u64(t11, 0x0F0F0F0F0F0F0F0F);
    const uint64_t t22 = _pdep_u64(t12, 0x0F0F0F0F0F0F0F0F);
    const __m128i t1 = _mm_set_epi64x(t22, t21);
    return _mm256_cvtepu8_epi16(t1);
}


static inline __m256i mirath_vector_extend_gf16_x32(const ff_t *in) {
    const uint16_t t11 = *((uint16_t *)(in + 0));
    const uint16_t t12 = *((uint16_t *)(in + 2));
    const uint16_t t13 = *((uint16_t *)(in + 4));
    const uint16_t t14 = *((uint16_t *)(in + 6));

    const uint64_t t21 = _pdep_u64(t11, 0x000F000F000F000F);
    const uint64_t t22 = _pdep_u64(t12, 0x000F000F000F000F);
    const uint64_t t23 = _pdep_u64(t13, 0x000F000F000F000F);
    const uint64_t t24 = _pdep_u64(t14, 0x000F000F000F000F);
    return _mm256_setr_epi64x(t21, t22, t23, t24);
}


#endif
