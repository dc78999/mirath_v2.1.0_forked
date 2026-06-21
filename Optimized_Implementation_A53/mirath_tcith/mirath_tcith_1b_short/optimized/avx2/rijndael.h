/**
 * @file rijndael.h
 * @brief Content for rijndael.h (AES-128 and Rijndael-256 optimized implementation)
 */

#ifndef RIJNDAEL_H
#define RIJNDAEL_H

#include <emmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <string.h>
#include <wmmintrin.h>

#define AES128_ROUNDS 10
#define RIJNDAEL256_ROUNDS 14

typedef struct {
    __m128i keys[AES128_ROUNDS + 1];
} aes_128_round_keys_t;

typedef struct {
    __m128i keys[(RIJNDAEL256_ROUNDS + 1)][2];
} rijndael_256_round_keys_t;

static inline __m128i aes_128_assist(__m128i temp1, __m128i temp2) {
    __m128i temp3;
    temp3 = _mm_slli_si128(temp1, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp1 = _mm_xor_si128(temp1, temp2);
    return temp1;
}

static inline void aes_128_key_expansion(aes_128_round_keys_t *round_keys, const unsigned char *key) {
    __m128i *Key_Schedule = round_keys->keys;
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);

    Key_Schedule[0] = *(__m128i *)key;

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    __m128i tmp = _mm_shuffle_epi8(Key_Schedule[0], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[1] = aes_128_assist(Key_Schedule[0], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[2] = aes_128_assist(Key_Schedule[1], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[2], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[3] = aes_128_assist(Key_Schedule[2], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[3], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[4] = aes_128_assist(Key_Schedule[3], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[4], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[5] = aes_128_assist(Key_Schedule[4], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[5], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[6] = aes_128_assist(Key_Schedule[5], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[6], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[7] = aes_128_assist(Key_Schedule[6], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[7], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[8] = aes_128_assist(Key_Schedule[7], tmp);

    rcon = _mm_set_epi32(0x1b, 0x1b, 0x1b, 0x1b);
    tmp = _mm_shuffle_epi8(Key_Schedule[8], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[9] = aes_128_assist(Key_Schedule[8], tmp);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[9], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    Key_Schedule[10] = aes_128_assist(Key_Schedule[9], tmp);
}

static inline void aes_128_encrypt(unsigned char *out, const unsigned char *in,
                                   const aes_128_round_keys_t *Key_Schedule) {
    __m128i const *KS = Key_Schedule->keys;
    __m128i data = _mm_load_si128(&((__m128i *)in)[0]);

    data = _mm_xor_si128(data, KS[0]);

    int j;
    for (j = 1; j < AES128_ROUNDS; j++) {
        data = _mm_aesenc_si128(data, KS[j]);
    }

    ((__m128i *)out)[0] = _mm_aesenclast_si128(data, KS[j]);
}

static inline void rijndael_256_assist(const __m128i round_key_in[2], __m128i temp1, __m128i round_key_out[2]) {
    __m128i t1, t2, t3, t4;

    t1 = round_key_in[0];
    t3 = round_key_in[1];
    t2 = temp1;

    t2 = _mm_shuffle_epi32(t2, 0xff);
    t4 = _mm_slli_si128(t1, 0x4);
    t1 = _mm_xor_si128(t1, t4);
    t4 = _mm_slli_si128(t4, 0x4);
    t1 = _mm_xor_si128(t1, t4);
    t4 = _mm_slli_si128(t4, 0x4);
    t1 = _mm_xor_si128(t1, t4);
    t1 = _mm_xor_si128(t1, t2);

    round_key_out[0] = t1;

    __m128i zero = {0};
    t4 = _mm_shuffle_epi32(t1, 0xFF);
    t4 = _mm_aesenclast_si128(t4, zero);

    t2 = _mm_shuffle_epi32(t4, 0xaa);
    t4 = _mm_slli_si128(t3, 0x4);
    t3 = _mm_xor_si128(t3, t4);
    t4 = _mm_slli_si128(t4, 0x4);
    t3 = _mm_xor_si128(t3, t4);
    t4 = _mm_slli_si128(t4, 0x4);
    t3 = _mm_xor_si128(t3, t4);
    t3 = _mm_xor_si128(t3, t2);

    round_key_out[1] = t3;
}

static inline void rijndael_256_key_expansion(rijndael_256_round_keys_t *round_keys, const unsigned char *key) {
    __m128i(*Key_Schedule)[2] = round_keys->keys;
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);

    Key_Schedule[0][0] = ((__m128i *)key)[0];
    Key_Schedule[0][1] = ((__m128i *)key)[1];

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    __m128i tmp = _mm_shuffle_epi8(Key_Schedule[0][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[0], tmp, Key_Schedule[1]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[1][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[1], tmp, Key_Schedule[2]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[2][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[2], tmp, Key_Schedule[3]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[3][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[3], tmp, Key_Schedule[4]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[4][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[4], tmp, Key_Schedule[5]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[5][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[5], tmp, Key_Schedule[6]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[6][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[6], tmp, Key_Schedule[7]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[7][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[7], tmp, Key_Schedule[8]);

    rcon = _mm_set_epi32(0x1B, 0x1B, 0x1B, 0x1B);
    tmp = _mm_shuffle_epi8(Key_Schedule[8][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[8], tmp, Key_Schedule[9]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[9][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[9], tmp, Key_Schedule[10]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[10][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[10], tmp, Key_Schedule[11]);

    rcon = _mm_slli_epi32(rcon, 1);
    tmp = _mm_shuffle_epi8(Key_Schedule[11][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[11], tmp, Key_Schedule[12]);

    rcon = _mm_set_epi32(0xAB, 0xAB, 0xAB, 0xAB);
    tmp = _mm_shuffle_epi8(Key_Schedule[12][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[12], tmp, Key_Schedule[13]);

    rcon = _mm_set_epi32(0x4D, 0x4D, 0x4D, 0x4D);
    tmp = _mm_shuffle_epi8(Key_Schedule[13][1], shuffle_mask);
    tmp = _mm_aesenclast_si128(tmp, rcon);
    rijndael_256_assist(Key_Schedule[13], tmp, Key_Schedule[14]);
}

static inline void rijndael_256_encrypt(unsigned char *out, const unsigned char *in, const rijndael_256_round_keys_t *Key_Schedule) {
    __m128i RIJNDAEL256_MASK = _mm_set_epi32(0x03020d0c, 0x0f0e0908, 0x0b0a0504, 0x07060100);
    __m128i BLEND_MASK = _mm_set_epi32(0x80000000, 0x80800000, 0x80800000, 0x80808000);

    __m128i const(*KS)[2] = Key_Schedule->keys;
    __m128i data1 = _mm_load_si128(&((__m128i *)in)[0]); /* load data block */
    __m128i data2 = _mm_load_si128(&((__m128i *)in)[1]);

    data1 = _mm_xor_si128(data1, KS[0][0]); /* round 0 (initial xor) */
    data2 = _mm_xor_si128(data2, KS[0][1]);

    int j;
    for (j = 1; j < RIJNDAEL256_ROUNDS; j++) {
        /*Blend to compensate for the shift rows shifts bytes between two
        128 bit blocks*/
        __m128i tmp1 = _mm_blendv_epi8(data1, data2, BLEND_MASK);
        __m128i tmp2 = _mm_blendv_epi8(data2, data1, BLEND_MASK);
        /*Shuffle that compensates for the additional shift in rows 3 and 4
        as opposed to rijndael128 (AES)*/
        tmp1 = _mm_shuffle_epi8(tmp1, RIJNDAEL256_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, RIJNDAEL256_MASK);
        /*This is the encryption step that includes sub bytes, shift rows,
        mix columns, xor with round key*/
        data1 = _mm_aesenc_si128(tmp1, KS[j][0]);
        data2 = _mm_aesenc_si128(tmp2, KS[j][1]);
    }
    __m128i tmp1 = _mm_blendv_epi8(data1, data2, BLEND_MASK);
    __m128i tmp2 = _mm_blendv_epi8(data2, data1, BLEND_MASK);
    tmp1 = _mm_shuffle_epi8(tmp1, RIJNDAEL256_MASK);
    tmp2 = _mm_shuffle_epi8(tmp2, RIJNDAEL256_MASK);
    ((__m128i *)out)[0] = _mm_aesenclast_si128(tmp1, KS[j][0]); /*last AES round */
    ((__m128i *)out)[1] = _mm_aesenclast_si128(tmp2, KS[j][1]);
}

#endif //RIJNDAEL_H
