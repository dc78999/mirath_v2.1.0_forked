#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <immintrin.h>

#include "data_type.h"

/// AES modulus
#define MODULUS 0x1B

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
}

static
const uint8_t __gf256_mulbase_avx2[256] __attribute__((aligned(32))) = {
        // repeated over each 128bit lane
        // 16**0 * i
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07, 0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f, 0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0,
        // 16**1 * i
        0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0, 0x00,0x1b,0x36,0x2d,0x6c,0x77,0x5a,0x41, 0xd8,0xc3,0xee,0xf5,0xb4,0xaf,0x82,0x99,
        // 16**2 * i
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e, 0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e, 0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb,
        // 16**3 * i
        0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb, 0x00,0x36,0x6c,0x5a,0xd8,0xee,0xb4,0x82, 0xab,0x9d,0xc7,0xf1,0x73,0x45,0x1f,0x29,
        // 16**4 * i
        0x00,0x04,0x08,0x0c,0x10,0x14,0x18,0x1c, 0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c, 0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed,
        // 16**5 * i
        0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed, 0x00,0x6c,0xd8,0xb4,0xab,0xc7,0x73,0x1f, 0x4d,0x21,0x95,0xf9,0xe6,0x8a,0x3e,0x52,
        // 16**6 * i
        0x00,0x08,0x10,0x18,0x20,0x28,0x30,0x38, 0x40,0x48,0x50,0x58,0x60,0x68,0x70,0x78, 0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1,
        // 16**7 * i
        0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1, 0x00,0xd8,0xab,0x73,0x4d,0x95,0xe6,0x3e, 0x9a,0x42,0x31,0xe9,0xd7,0x0f,0x7c,0xa4
};

//// extended for avx2 registers
static
const uint8_t __gf256_mulbase_avx_v2[256] __attribute__((aligned(32))) = {
        // 1*i
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07, 0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f, 0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0,
        // 2*i
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e, 0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e, 0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb,
        // 4*i
        0x00,0x04,0x08,0x0c,0x10,0x14,0x18,0x1c, 0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c, 0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed,
        // 8*i
        0x00,0x08,0x10,0x18,0x20,0x28,0x30,0x38, 0x40,0x48,0x50,0x58,0x60,0x68,0x70,0x78, 0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1,
        // 16*i
        0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0, 0x00,0x1b,0x36,0x2d,0x6c,0x77,0x5a,0x41, 0xd8,0xc3,0xee,0xf5,0xb4,0xaf,0x82,0x99,
        // 32*i
        0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb, 0x00,0x36,0x6c,0x5a,0xd8,0xee,0xb4,0x82, 0xab,0x9d,0xc7,0xf1,0x73,0x45,0x1f,0x29,
        // 64*i
        0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed, 0x00,0x6c,0xd8,0xb4,0xab,0xc7,0x73,0x1f, 0x4d,0x21,0x95,0xf9,0xe6,0x8a,0x3e,0x52,
        // 128*i
        0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1, 0x00,0xd8,0xab,0x73,0x4d,0x95,0xe6,0x3e, 0x9a,0x42,0x31,0xe9,0xd7,0x0f,0x7c,0xa4
};


/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    ff_mu_t result = -(a & 1) & b;
    ff_mu_t tmp = b;
    for(int i=1 ; i<8 ; i++) {
        tmp = (tmp << 1) ^ (-(tmp >> 7) & MODULUS);
        result = result ^ (-(a >> i & 1) & tmp);
    }
    return result;
}

/// NOTE: is actually constant time
/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    ff_mu_t result = a;
    for(int i=0 ; i<6 ; i++) {
        result = mirath_ff_mu_mult(result, result);
        result = mirath_ff_mu_mult(result, a);
    }
    result = mirath_ff_mu_mult(result, result);
    return result;
}

/// horizontal add, but not withing a single limb, but over the 32x 8-bit limbs
/// NOTE: add = xor
/// \param in[in]: [in_0, ..., in_31]
/// \return sum(in_0, in_1, ..., in[32])
static inline uint32_t mirath_ff_mu_hadd_u32_u256(const __m256i in) {
    __m256i ret = _mm256_xor_si256(in, _mm256_srli_si256(in, 4));
    ret = _mm256_xor_si256(ret, _mm256_srli_si256(ret, 8));
    ret = _mm256_xor_si256(ret, _mm256_permute2x128_si256(ret, ret, 129)); // 0b10000001
    return _mm256_extract_epi32(ret, 0);
}

/// \param a = [a_0, ..., a_31]
/// \param b = [b_0, ..., b_31]
/// \return [a_0*b_0, ..., a_31*b_31]
static inline
__m256i gf256v_mul_u256(const __m256i a,
                        const __m256i b) {
    const __m256i zero = _mm256_set1_epi32(0),
            mask = _mm256_set1_epi8(0x1B);
    __m256i r;

    r = _mm256_blendv_epi8(zero, a, b);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 1)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 2)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 3)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 4)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 5)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 6)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    r = _mm256_blendv_epi8(zero, a, _mm256_slli_epi16(b, 7)) ^ _mm256_blendv_epi8(zero, mask, r) ^ _mm256_add_epi8(r, r);
    return r;
}

/// SRC: https://github.com/pqov/pqov-paper
/// NOTE: this is actually just a scalar multiplication, as `tab_l` and `tab_h`
///     are computed from a single scalar.
/// \param tab_l[in]: precomputation from `mirath_ff_mu_generate_multab_16_single_element_u256`
/// \param tab_h[in]: precomputation from `mirath_ff_mu_generate_multab_16_single_element_u256`
/// \param v[in]: [tab_l, tab_h] * v
/// \param mask_f[in] = [0xF, 0xF, ..., 0xF]
/// \return [tab_l, tab_h] * v
static inline
__m256i mirath_ff_mu_linear_transform_8x8_256b(__m256i tab_l,
                                        __m256i tab_h,
                                        __m256i v,
                                        __m256i mask_f) {
    return _mm256_shuffle_epi8(tab_l, v & mask_f) ^
           _mm256_shuffle_epi8(tab_h, _mm256_srli_epi16(v, 4) & mask_f);
}

/// SRC: https://github.com/pqov/pqov-paper
/// \param tab_l[in]: precomputation from `mirath_ff_mu_generate_multab_16_single_element_u256`
/// \param tab_h[in]: precomputation from `mirath_ff_mu_generate_multab_16_single_element_u256`
/// \param v[in]: [tab_l, tab_h] * v
/// \param mask_f[in] = [0xF, 0xF, ..., 0xF]
/// \return [tab_l, tab_h] * v
static inline
__m128i mirath_ff_mu_linear_transform_8x8_128b(__m128i tab_l,
                                        __m128i tab_h,
                                        __m128i v,
                                        __m128i mask_f) {
    return _mm_shuffle_epi8(tab_l, v & mask_f) ^
           _mm_shuffle_epi8(tab_h, _mm_srli_epi16(v, 4) & mask_f);
}

/// SRC: https://github.com/pqov/pqov-paper
/// \param a[in]: scalar
/// \return a precomputation table for scalar*v, where v is a avx register
static inline
__m256i mirath_ff_mu_generate_multab_16_single_element_u256(const uint8_t a) {
    __m256i bx = _mm256_set1_epi16(a);
    __m256i b1 = _mm256_srli_epi16(bx , 1 );

    __m256i tab0 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*0));
    __m256i tab1 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*1));
    __m256i tab2 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*2));
    __m256i tab3 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*3));
    __m256i tab4 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*4));
    __m256i tab5 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*5));
    __m256i tab6 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*6));
    __m256i tab7 = _mm256_load_si256((__m256i const *) (__gf256_mulbase_avx_v2 + 32*7));

    __m256i mask_1  = _mm256_set1_epi16(1);
    __m256i mask_4  = _mm256_set1_epi16(4);
    __m256i mask_16 = _mm256_set1_epi16(16);
    __m256i mask_64 = _mm256_set1_epi16(64);
    __m256i mask_0  = _mm256_setzero_si256();

    return ( tab0 & _mm256_cmpgt_epi16( bx&mask_1  , mask_0) )
           ^ ( tab1 & _mm256_cmpgt_epi16( b1&mask_1  , mask_0) )
           ^ ( tab2 & _mm256_cmpgt_epi16( bx&mask_4  , mask_0) )
           ^ ( tab3 & _mm256_cmpgt_epi16( b1&mask_4  , mask_0) )
           ^ ( tab4 & _mm256_cmpgt_epi16( bx&mask_16 , mask_0) )
           ^ ( tab5 & _mm256_cmpgt_epi16( b1&mask_16 , mask_0) )
           ^ ( tab6 & _mm256_cmpgt_epi16( bx&mask_64 , mask_0) )
           ^ ( tab7 & _mm256_cmpgt_epi16( b1&mask_64 , mask_0) );
}

/// \param in[in] 4 bytes
///     in = [in_0, ..., in_7, in_8, ..., in_31] in bits
/// \return [in_0, ..., in_31], extended to 8-bit limbs
static inline __m256i mirath_expand_gf2_x32_u256(const uint8_t *in) {
    const uint8_t t11 = *(in + 0);
    const uint8_t t12 = *(in + 1);
    const uint8_t t13 = *(in + 2);
    const uint8_t t14 = *(in + 3);

    const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
    const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
    const uint64_t t23 = _pdep_u64(t13, 0x0101010101010101);
    const uint64_t t24 = _pdep_u64(t14, 0x0101010101010101);
    return _mm256_setr_epi64x(t21, t22, t23, t24);
}

/// \param in[in] 2 bytes
///     in = [in_0, ..., in_7, in_8, ..., in_15] in bitsk
/// \return [in_0, ..., in_15], extended to 8-bit limbs
static inline __m128i mirath_expand_gf2_x16_u256(const uint8_t *in) {
    const uint32_t t11 = *(in + 0);
    const uint32_t t12 = *(in + 1);
    const uint64_t t21 = _pdep_u64(t11, 0x0101010101010101);
    const uint64_t t22 = _pdep_u64(t12, 0x0101010101010101);
    return _mm_set_epi64x(t22, t21);
}

/// \param in[in] 1 bytes
///     in = [in_0, ..., in_7] in bitsk
/// \return [in_0, ..., in_7], extended to 8-bit limbs
static inline uint64_t mirath_expand_gf2_x8_u256(const uint8_t *in) {
    const uint32_t t11 = *(in + 0);
    return _pdep_u64(t11, 0x0101010101010101);
}
/// \param a in gf2to8
/// \param b in gf2, not compresses: a single bit in
/// \return a*b
static inline __m256i mirath_ff_mu_mul_ff_u256(const __m256i a,
                                               const __m256i b) {
    const __m256i m1 = _mm256_set1_epi8(-1);
    const __m256i t1 = _mm256_sign_epi8(b, m1);
    return a & t1;
}

/// \param a in gf2to8
/// \param b in gf2, not compresses: a single bit in
/// \return a*b
static inline __m128i mirath_ff_mu_mul_ff_u128(const __m128i a,
                                             const __m128i b) {
    const __m128i m1 = _mm_set1_epi8(-1);
    const __m128i t1 = _mm_sign_epi8(b, m1);
    return a & t1;
}

#endif
