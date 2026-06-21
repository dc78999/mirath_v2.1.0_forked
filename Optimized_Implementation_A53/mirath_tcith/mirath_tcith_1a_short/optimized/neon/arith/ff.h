/**
 * @file ff.h
 * @brief GF(2^4) finite field arithmetic with ARM NEON optimizations.
 *
 * Constant-time implementation: all operations use fixed-time NEON instructions.
 * No secret-dependent branches or memory accesses.
 *
 * NEON intrinsics used:
 *   vqtbl1q_u8  - register-based table lookup (constant-time, no cache access)
 *   vbslq_u8    - bitwise select (constant-time conditional move)
 *   veorq_u8    - XOR
 *   vshlq_n_u16 - shift left
 */

#ifndef ARITH_FF_H
#define ARITH_FF_H

#include <stdint.h>
#include <arm_neon.h>

#include "data_type.h"

#define FF_MODULUS 3u

/// i*j = [i*16 + j] — full multiplication table for GF(2^4)
static const ff_t mirath_ff_mult_table[256] = {
    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
    0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
    0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e,0x03,0x01,0x07,0x05,0x0b,0x09,0x0f,0x0d,
    0x00,0x03,0x06,0x05,0x0c,0x0f,0x0a,0x09,0x0b,0x08,0x0d,0x0e,0x07,0x04,0x01,0x02,
    0x00,0x04,0x08,0x0c,0x03,0x07,0x0b,0x0f,0x06,0x02,0x0e,0x0a,0x05,0x01,0x0d,0x09,
    0x00,0x05,0x0a,0x0f,0x07,0x02,0x0d,0x08,0x0e,0x0b,0x04,0x01,0x09,0x0c,0x03,0x06,
    0x00,0x06,0x0c,0x0a,0x0b,0x0d,0x07,0x01,0x05,0x03,0x09,0x0f,0x0e,0x08,0x02,0x04,
    0x00,0x07,0x0e,0x09,0x0f,0x08,0x01,0x06,0x0d,0x0a,0x03,0x04,0x02,0x05,0x0c,0x0b,
    0x00,0x08,0x03,0x0b,0x06,0x0e,0x05,0x0d,0x0c,0x04,0x0f,0x07,0x0a,0x02,0x09,0x01,
    0x00,0x09,0x01,0x08,0x02,0x0b,0x03,0x0a,0x04,0x0d,0x05,0x0c,0x06,0x0f,0x07,0x0e,
    0x00,0x0a,0x07,0x0d,0x0e,0x04,0x09,0x03,0x0f,0x05,0x08,0x02,0x01,0x0b,0x06,0x0c,
    0x00,0x0b,0x05,0x0e,0x0a,0x01,0x0f,0x04,0x07,0x0c,0x02,0x09,0x0d,0x06,0x08,0x03,
    0x00,0x0c,0x0b,0x07,0x05,0x09,0x0e,0x02,0x0a,0x06,0x01,0x0d,0x0f,0x03,0x04,0x08,
    0x00,0x0d,0x09,0x04,0x01,0x0c,0x08,0x05,0x02,0x0f,0x0b,0x06,0x03,0x0e,0x0a,0x07,
    0x00,0x0e,0x0f,0x01,0x0d,0x03,0x02,0x0c,0x09,0x07,0x06,0x08,0x04,0x0a,0x0b,0x05,
    0x00,0x0f,0x0d,0x02,0x09,0x06,0x04,0x0b,0x01,0x0e,0x0c,0x03,0x08,0x07,0x05,0x0a,
};

/// M[i, j] = 16**i * j mod 16 for i = row, j = column
/// Aligned for NEON loads. Each row is 16 bytes (one uint8x16_t).
static const uint8_t mirath_ff_mulbase[64] __attribute__((aligned(16))) = {
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07, 0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e, 0x03,0x01,0x07,0x05,0x0b,0x09,0x0f,0x0d,
        0x00,0x04,0x08,0x0c,0x03,0x07,0x0b,0x0f, 0x06,0x02,0x0e,0x0a,0x05,0x01,0x0d,0x09,
        0x00,0x08,0x03,0x0b,0x06,0x0e,0x05,0x0d, 0x0c,0x04,0x0f,0x07,0x0a,0x02,0x09,0x01,
};

/// i*{-1} = [i]
static const ff_t mirath_ff_inv_table[16] __attribute__((aligned(16))) = {
    0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8
};

/// \return a+b mod q  (XOR in GF(2^4))
static inline ff_t mirath_ff_add(const ff_t a, const ff_t b) {
    return a ^ b;
}

/// NOTE: assumes a mod 16
/// \return a^{-1} — only use on PUBLIC indices for constant-time safety
static inline ff_t mirath_ff_inv(const ff_t a) {
    return mirath_ff_inv_table[a];
}

/// Constant-time GF(2^4) multiplication using bitwise operations.
/// No branches, no secret-dependent memory accesses.
/// \return a*b mod (x^4 + x + 1)
static inline ff_t mirath_ff_product(const ff_t a, const ff_t b) {
    uint8_t r;
    r =    (-(b>>3    ) & a);
    r =    (-(b>>2 & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
    r =    (-(b>>1 & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
    return (-(b    & 1) & a) ^ (-(r>>3) & FF_MODULUS) ^ ((r+r) & 0x0F);
}

/// Horizontal XOR reduction over a 128-bit NEON register, returning a uint64_t.
/// Constant-time: all NEON integer operations on A53 are fixed-latency.
static inline uint64_t mirath_hadd_u128_64(const uint8x16_t in) {
    /* XOR the upper and lower 64-bit halves */
    uint64x2_t v = vreinterpretq_u64_u8(in);
    return vgetq_lane_u64(v, 0) ^ vgetq_lane_u64(v, 1);
}

/// NEON GF(2^4) multiplication: multiply 32 nibbles (16 bytes) at once.
/// Uses vqtbl1q_u8 for register-based lookup (constant-time, no cache).
///
/// Each byte contains two nibbles (low and high). We process them separately
/// using table lookups for the multiplication-by-power tables.
///
/// \param a packed nibble vector (16 bytes = 32 nibbles)
/// \param _b packed nibble vector
/// \return a*b element-wise in GF(2^4), packed
static inline uint8x16_t mirath_ff_mul_full_u128(const uint8x16_t a,
                                                  const uint8x16_t _b) {
    const uint8x16_t mask_lvl2 = vld1q_u8(mirath_ff_mulbase + 16);
    const uint8x16_t mask_lvl3 = vld1q_u8(mirath_ff_mulbase + 32);
    const uint8x16_t mask_lvl4 = vld1q_u8(mirath_ff_mulbase + 48);
    const uint8x16_t zero = vdupq_n_u8(0);
    const uint8x16_t mask1 = vdupq_n_u8(0x0F);

    const uint8x16_t b = vandq_u8(_b, mask1);
    const uint8x16_t b2 = vandq_u8(vshrq_n_u8(vreinterpretq_u8_u16(vshrq_n_u16(vreinterpretq_u16_u8(_b), 4)), 0), mask1);

    /* Level 0: identity table */
    uint8x16_t low_lookup  = b;
    uint8x16_t high_lookup = vshlq_n_u8(b2, 4);

    /* Bit 3 of a selects whether to include the lookup result */
    uint8x16_t sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 7 - 3)), 7));
    uint8x16_t sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 3 - 3)), 7));
    uint8x16_t tmp = veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                               vbslq_u8(sel_hi, high_lookup, zero));

    /* Level 1 */
    low_lookup = vqtbl1q_u8(mask_lvl2, b);
    high_lookup = vshlq_n_u8(vqtbl1q_u8(mask_lvl2, b2), 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 7 - 2)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 3 - 2)), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    /* Level 2 */
    low_lookup = vqtbl1q_u8(mask_lvl3, b);
    high_lookup = vshlq_n_u8(vqtbl1q_u8(mask_lvl3, b2), 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 7 - 1)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 3 - 1)), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    /* Level 3 */
    low_lookup = vqtbl1q_u8(mask_lvl4, b);
    high_lookup = vshlq_n_u8(vqtbl1q_u8(mask_lvl4, b2), 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 7 - 0)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(a), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    return tmp;
}

/// NEON GF(2^4) multiplication where both nibbles in each byte of b are identical.
/// Simplified version of mirath_ff_mul_full_u128.
static inline uint8x16_t mirath_ff_mul_u128(const uint8x16_t a,
                                             const uint8x16_t b) {
    const uint8x16_t mask_lvl2 = vld1q_u8(mirath_ff_mulbase + 16);
    const uint8x16_t mask_lvl3 = vld1q_u8(mirath_ff_mulbase + 32);
    const uint8x16_t mask_lvl4 = vld1q_u8(mirath_ff_mulbase + 48);
    const uint8x16_t zero = vdupq_n_u8(0);

    uint8x16_t low_lookup = b;
    uint8x16_t high_lookup = vshlq_n_u8(low_lookup, 4);

    /* Bit 3 of a */
    uint8x16_t sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 4)), 7));
    uint8x16_t sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(a), 7));
    uint8x16_t tmp = veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                               vbslq_u8(sel_hi, high_lookup, zero));

    /* Bit 2 */
    low_lookup = vqtbl1q_u8(mask_lvl2, b);
    high_lookup = vshlq_n_u8(low_lookup, 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 5)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 1)), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    /* Bit 1 */
    low_lookup = vqtbl1q_u8(mask_lvl3, b);
    high_lookup = vshlq_n_u8(low_lookup, 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 6)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 2)), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    /* Bit 0 */
    low_lookup = vqtbl1q_u8(mask_lvl4, b);
    high_lookup = vshlq_n_u8(low_lookup, 4);
    sel_lo = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 7)), 7));
    sel_hi = vreinterpretq_u8_s8(vshrq_n_s8(vreinterpretq_s8_u8(vshlq_n_u8(a, 3)), 7));
    tmp = veorq_u8(tmp, veorq_u8(vbslq_u8(sel_lo, low_lookup, zero),
                                  vbslq_u8(sel_hi, high_lookup, zero)));

    return tmp;
}

/// Generate multiplication table for a single GF(2^4) element (NEON 128-bit).
/// Result: a uint8x16_t where entry i = b * i for i in 0..15.
/// Uses constant-time comparison + bitwise select.
static inline
uint8x16_t mirath_ff_tbl16_multab(const ff_t b) {
    uint16x8_t bx = vdupq_n_u16(b & 0xf);
    uint16x8_t b1 = vshrq_n_u16(bx, 1);

    uint8x16_t tab0 = vld1q_u8(mirath_ff_mulbase + 16*0);
    uint8x16_t tab1 = vld1q_u8(mirath_ff_mulbase + 16*1);
    uint8x16_t tab2 = vld1q_u8(mirath_ff_mulbase + 16*2);
    uint8x16_t tab3 = vld1q_u8(mirath_ff_mulbase + 16*3);

    uint16x8_t mask_1 = vdupq_n_u16(1);
    uint16x8_t mask_4 = vdupq_n_u16(4);
    uint16x8_t mask_0 = vdupq_n_u16(0);

    /* Constant-time: vcgtq_u16 + vandq produce masks without branches */
    uint16x8_t s0 = vcgtq_u16(vandq_u16(bx, mask_1), mask_0);
    uint16x8_t s1 = vcgtq_u16(vandq_u16(b1, mask_1), mask_0);
    uint16x8_t s2 = vcgtq_u16(vandq_u16(bx, mask_4), mask_0);
    uint16x8_t s3 = vcgtq_u16(vandq_u16(b1, mask_4), mask_0);

    uint8x16_t result = vandq_u8(tab0, vreinterpretq_u8_u16(s0));
    result = veorq_u8(result, vandq_u8(tab1, vreinterpretq_u8_u16(s1)));
    result = veorq_u8(result, vandq_u8(tab2, vreinterpretq_u8_u16(s2)));
    result = veorq_u8(result, vandq_u8(tab3, vreinterpretq_u8_u16(s3)));

    return result;
}

/// NEON linear transform: lookup table multiplication.
/// Equivalent to _mm_shuffle_epi8 pattern used in AVX2.
/// Constant-time: vqtbl1q_u8 is register-based.
static inline
uint8x16_t mirath_ff_linear_transform_8x8_128b(uint8x16_t tab_l,
                                                 uint8x16_t tab_h,
                                                 uint8x16_t v,
                                                 uint8x16_t mask_f) {
    return veorq_u8(vqtbl1q_u8(tab_l, vandq_u8(v, mask_f)),
                    vqtbl1q_u8(tab_h, vandq_u8(vshrq_n_u8(v, 4), mask_f)));
}

#endif /* ARITH_FF_H */
