#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <immintrin.h>

#include "data_type.h"

/// AES modulus
#define MODULUS 0x1B
#define MASK_LSB_PER_BIT ((uint64_t)0x0101010101010101)
#define MASK_MSB_PER_BIT (MASK_LSB_PER_BIT*0x80)
#define MASK_XLSB_PER_BIT (MASK_LSB_PER_BIT*0xFE)

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

static
const uint8_t __gf256_mulbase_avx_v2[256] __attribute__((aligned(32))) = {
        // 1 * i
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07, 0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f, 0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0,
        // 2 * i
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e, 0x10,0x12,0x14,0x16,0x18,0x1a,0x1c,0x1e, 0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb,
        // 4 * i
        0x00,0x04,0x08,0x0c,0x10,0x14,0x18,0x1c, 0x20,0x24,0x28,0x2c,0x30,0x34,0x38,0x3c, 0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed,
        // 8 * i
        0x00,0x08,0x10,0x18,0x20,0x28,0x30,0x38, 0x40,0x48,0x50,0x58,0x60,0x68,0x70,0x78, 0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1,
        // 16 * i
        0x00,0x10,0x20,0x30,0x40,0x50,0x60,0x70, 0x80,0x90,0xa0,0xb0,0xc0,0xd0,0xe0,0xf0, 0x00,0x1b,0x36,0x2d,0x6c,0x77,0x5a,0x41, 0xd8,0xc3,0xee,0xf5,0xb4,0xaf,0x82,0x99,
        // 32 * i
        0x00,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, 0x1b,0x3b,0x5b,0x7b,0x9b,0xbb,0xdb,0xfb, 0x00,0x36,0x6c,0x5a,0xd8,0xee,0xb4,0x82, 0xab,0x9d,0xc7,0xf1,0x73,0x45,0x1f,0x29,
        // 64 * i
        0x00,0x40,0x80,0xc0,0x1b,0x5b,0x9b,0xdb, 0x36,0x76,0xb6,0xf6,0x2d,0x6d,0xad,0xed, 0x00,0x6c,0xd8,0xb4,0xab,0xc7,0x73,0x1f, 0x4d,0x21,0x95,0xf9,0xe6,0x8a,0x3e,0x52,
        // 128 * i
        0x00,0x80,0x1b,0x9b,0x36,0xb6,0x2d,0xad, 0x6c,0xec,0x77,0xf7,0x5a,0xda,0x41,0xc1, 0x00,0xd8,0xab,0x73,0x4d,0x95,0xe6,0x3e, 0x9a,0x42,0x31,0xe9,0xd7,0x0f,0x7c,0xa4
};

// size doubled for avx, so its of size 32, which is one register
static const uint8_t mirath_map_ff_to_ff_mu[32] __attribute__((aligned(32))) = {
        0, 1, 92, 93, 224, 225, 188, 189, 80, 81, 12, 13, 176, 177, 236, 237,
        0, 1, 92, 93, 224, 225, 188, 189, 80, 81, 12, 13, 176, 177, 236, 237
};

// Warning, getting the inverse of a secret value using this table,
//   would lead to non-constant-time implementation. Only accessing
//   to public positions is allowed.
static const uint8_t mirath_ff_mu_inv_table[256] = {0, 1, 141, 246, 203, 82, 123, 209, 232, 79, 41, 192, 176, 225, 229, 199, 116, 180, 170, 75, 153, 43, 96, 95, 88, 63, 253, 204, 255, 64, 238, 178, 58, 110, 90, 241, 85, 77, 168, 201, 193, 10, 152, 21, 48, 68, 162, 194, 44, 69, 146, 108, 243, 57, 102, 66, 242, 53, 32, 111, 119, 187, 89, 25, 29, 254, 55, 103, 45, 49, 245, 105, 167, 100, 171, 19, 84, 37, 233, 9, 237, 92, 5, 202, 76, 36, 135, 191, 24, 62, 34, 240, 81, 236, 97, 23, 22, 94, 175, 211, 73, 166, 54, 67, 244, 71, 145, 223, 51, 147, 33, 59, 121, 183, 151, 133, 16, 181, 186, 60, 182, 112, 208, 6, 161, 250, 129, 130, 131, 126, 127, 128, 150, 115, 190, 86, 155, 158, 149, 217, 247, 2, 185, 164, 222, 106, 50, 109, 216, 138, 132, 114, 42, 20, 159, 136, 249, 220, 137, 154, 251, 124, 46, 195, 143, 184, 101, 72, 38, 200, 18, 74, 206, 231, 210, 98, 12, 224, 31, 239, 17, 117, 120, 113, 165, 142, 118, 61, 189, 188, 134, 87, 11, 40, 47, 163, 218, 212, 228, 15, 169, 39, 83, 4, 27, 252, 172, 230, 122, 7, 174, 99, 197, 219, 226, 234, 148, 139, 196, 213, 157, 248, 144, 107, 177, 13, 214, 235, 198, 14, 207, 173, 8, 78, 215, 227, 93, 80, 30, 179, 91, 35, 56, 52, 104, 70, 3, 140, 221, 156, 125, 160, 205, 26, 65, 28};

/// \return a+b
static inline ff_mu_t mirath_ff_mu_add(const ff_mu_t a, const ff_mu_t b) {
    return a^b;
}

/// \return a*b
static inline ff_mu_t mirath_ff_mu_mult(const ff_mu_t a, const ff_mu_t b) {
    ff_mu_t r;
    r = (-(b>>7u     ) & a);
    r = (-(b>>6u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    r = (-(b>>5u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    r = (-(b>>4u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    r = (-(b>>3u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    r = (-(b>>2u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    r = (-(b>>1u & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
    return (-(b  & 1u) & a) ^ (-(r>>7u) & MODULUS) ^ (r+r);
}

/// \return a^-1
static inline ff_mu_t mirath_ff_mu_inv(const ff_mu_t a) {
    return mirath_ff_mu_inv_table[a];
}

/// loads 32 elements/16 bytes from gf16 and extends them to gf256
static inline __m256i mirath_ff_mu_extend_gf16_x32(const uint8_t *in) {
    const __m128i tbl   = _mm_load_si128((__m128i *)mirath_map_ff_to_ff_mu);
    const __m128i mask  = _mm_set1_epi8(0x0F);

    const __m128i load = _mm_loadu_si128((__m128i *) (in));
    const __m128i ll = _mm_shuffle_epi8(tbl, load & mask);
    const __m128i lh = _mm_shuffle_epi8(tbl, _mm_srli_epi16(load, 4) & mask);
    const __m256i tl = _mm256_setr_m128i(ll, _mm_bsrli_si128(ll, 8));
    const __m256i th = _mm256_setr_m128i(lh, _mm_bsrli_si128(lh, 8));
    return _mm256_unpacklo_epi8(tl, th);
}


/// horizontal xor, but not withing a single limb, but over the 8 -32bit limbs
/// NOTE: adding is just a xor.
/// \param in[in] avx register with 8 32-bit limbs. 
/// \return sum [in_0, in_1, ..., in_8]
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
#endif
