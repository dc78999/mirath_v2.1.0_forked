#ifndef ARITH_FF_MU_H
#define ARITH_FF_MU_H

#include <stdint.h>
#include <arm_neon.h>

#include "data_type.h"

/// AES modulus
#define MODULUS 0x1B
#define MASK_LSB_PER_BIT ((uint64_t)0x0101010101010101)
#define MASK_MSB_PER_BIT (MASK_LSB_PER_BIT*0x80)
#define MASK_XLSB_PER_BIT (MASK_LSB_PER_BIT*0xFE)

static inline void mirath_vector_set_to_ff_mu(ff_mu_t *v, const uint8_t *sample, const uint32_t n) {
    memcpy(v, sample, n);
}

// 256 bytes total. For each row: lower 16 bytes is tab_l, upper 16 bytes is tab_h
static const uint8_t __gf256_mulbase_neon[256] __attribute__((aligned(16))) = {
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


static const uint8_t mirath_map_ff_to_ff_mu[16] __attribute__((aligned(16))) = {
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

/// loads 16 elements/8 bytes from gf16 and extends them to gf256
static inline uint8x16_t mirath_ff_mu_extend_gf16_x16(const uint8_t *in) {
    uint8x8_t load8 = vld1_u8(in);
    uint8x16_t load16 = vcombine_u8(load8, vdup_n_u8(0));
    uint8x16_t mask = vdupq_n_u8(0x0F);
    uint8x16_t tbl = vld1q_u8(mirath_map_ff_to_ff_mu);
    
    uint8x16_t ll = vqtbl1q_u8(tbl, vandq_u8(load16, mask));
    uint8x16_t lh = vqtbl1q_u8(tbl, vandq_u8(vshrq_n_u8(load16, 4), mask));
    
    uint8x8x2_t res = vzip_u8(vget_low_u8(ll), vget_low_u8(lh));
    return vcombine_u8(res.val[0], res.val[1]);
}

/// horizontal xor, over 4 32bit limbs
static inline uint32_t mirath_ff_mu_hadd_u32_u128(const uint8x16_t in) {
    uint32x4_t v32 = vreinterpretq_u32_u8(in);
    uint32x2_t low = vget_low_u32(v32);
    uint32x2_t high = vget_high_u32(v32);
    uint32x2_t sum = veor_u32(low, high);
    sum = vpadd_u32(sum, sum);
    return vget_lane_u32(sum, 0);
}

static inline
uint8x16_t gf256v_mul_u128(const uint8x16_t a,
                           const uint8x16_t b) {
    const uint8x16_t zero = vdupq_n_u8(0);
    const uint8x16_t mask = vdupq_n_u8(0x1B);
    uint8x16_t r;

    r = vbslq_u8(vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b), 7)), a, zero);
    
    uint8x16_t b1 = vshlq_n_u8(b, 1);
    uint8x16_t mask1 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b1), 7));
    uint8x16_t msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask1, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b2 = vshlq_n_u8(b, 2);
    uint8x16_t mask2 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b2), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask2, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b3 = vshlq_n_u8(b, 3);
    uint8x16_t mask3 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b3), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask3, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b4 = vshlq_n_u8(b, 4);
    uint8x16_t mask4 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b4), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask4, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b5 = vshlq_n_u8(b, 5);
    uint8x16_t mask5 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b5), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask5, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b6 = vshlq_n_u8(b, 6);
    uint8x16_t mask6 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b6), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask6, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    uint8x16_t b7 = vshlq_n_u8(b, 7);
    uint8x16_t mask7 = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(b7), 7));
    msb_r = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(r), 7));
    r = veorq_u8(veorq_u8(vbslq_u8(mask7, a, zero), vbslq_u8(msb_r, mask, zero)), vshlq_n_u8(r, 1));
    
    return r;
}

static inline
uint8x16_t mirath_ff_mu_linear_transform_8x8_128b(uint8x16_t tab_l,
                                                 uint8x16_t tab_h,
                                                 uint8x16_t v,
                                                 uint8x16_t mask_f) {
    return veorq_u8(vqtbl1q_u8(tab_l, vandq_u8(v, mask_f)),
                    vqtbl1q_u8(tab_h, vandq_u8(vshrq_n_u8(v, 4), mask_f)));
}

// Returns a pair of uint8x16_t containing tab_l and tab_h
static inline
uint8x16x2_t mirath_ff_mu_generate_multab_16_single_element_u128(const uint8_t a) {
    uint16x8_t bx = vdupq_n_u16(a);
    uint16x8_t b1 = vshrq_n_u16(bx , 1);

    uint16x8_t mask_1  = vdupq_n_u16(1);
    uint16x8_t mask_4  = vdupq_n_u16(4);
    uint16x8_t mask_16 = vdupq_n_u16(16);
    uint16x8_t mask_64 = vdupq_n_u16(64);
    uint16x8_t mask_0  = vdupq_n_u16(0);

    uint8x16_t m0 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(bx, mask_1), mask_0));
    uint8x16_t m1 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(b1, mask_1), mask_0));
    uint8x16_t m2 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(bx, mask_4), mask_0));
    uint8x16_t m3 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(b1, mask_4), mask_0));
    uint8x16_t m4 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(bx, mask_16), mask_0));
    uint8x16_t m5 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(b1, mask_16), mask_0));
    uint8x16_t m6 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(bx, mask_64), mask_0));
    uint8x16_t m7 = vreinterpretq_u8_u16(vcgtq_u16(vandq_u16(b1, mask_64), mask_0));

    uint8x16x2_t res;

    // generate tab_l
    res.val[0] = vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*0), m0);
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*1), m1));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*2), m2));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*3), m3));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*4), m4));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*5), m5));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*6), m6));
    res.val[0] = veorq_u8(res.val[0], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*7), m7));

    // generate tab_h
    res.val[1] = vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*0 + 16), m0);
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*1 + 16), m1));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*2 + 16), m2));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*3 + 16), m3));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*4 + 16), m4));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*5 + 16), m5));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*6 + 16), m6));
    res.val[1] = veorq_u8(res.val[1], vandq_u8(vld1q_u8(__gf256_mulbase_neon + 32*7 + 16), m7));

    return res;
}

#endif
