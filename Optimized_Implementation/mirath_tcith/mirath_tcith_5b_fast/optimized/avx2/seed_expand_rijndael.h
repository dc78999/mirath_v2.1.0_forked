/**
 * @file seed_expand_functions_avx.h
 * @brief Content for seed_expand_functions_avx.h (Seed expand functions based on Rijndael-256)
 */

#ifndef MIRATH_SEED_EXPAND_RIJNDAEL_256_H
#define MIRATH_SEED_EXPAND_RIJNDAEL_256_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

typedef uint8_t block256_t[32] __attribute__ ((aligned (32)));

static inline void rijndael_expand_seed(uint8_t dst[2][32], const uint8_t salt[32], const uint32_t idx, const uint8_t seed[32]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    __m128i block_0[2] = {0};
    __m128i block_1[2] = {0};

    block256_t *const output = dst;
    const uint8_t *const seed_with_zeros = seed;

    uint8_t *msg = (uint8_t *)block_0;

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 32);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;

    // salt ^ (domain_separator || idx || 1)
    block_1[0] = block_0[0];
    block_1[1] = block_0[1];
    msg = (uint8_t *)block_1;
    msg[0] ^= 0x01;

    __m128i round_key[2];
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);
    __m128i RIJNDAEL256_MASK = _mm_set_epi32(0x03020d0c, 0x0f0e0908, 0x0b0a0504, 0x07060100);
    __m128i BLEND_MASK = _mm_set_epi32(0x80000000, 0x80800000, 0x80800000, 0x80808000);

    round_key[0] = ((__m128i *)seed_with_zeros)[0];
    round_key[1] = ((__m128i *)seed_with_zeros)[1];

    block_0[0] = _mm_xor_si128(block_0[0], round_key[0]); /* round 0 (initial xor) */
    block_0[1] = _mm_xor_si128(block_0[1], round_key[1]);
    block_1[0] = _mm_xor_si128(block_1[0], round_key[0]); /* round 0 (initial xor) */
    block_1[1] = _mm_xor_si128(block_1[1], round_key[1]);

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    for (int i = 1; i < 9; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);
        rcon = _mm_slli_epi32(rcon, 1);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0x1B, 0x1B, 0x1B, 0x1B);
    for (int i = 9; i < 13; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);
        rcon = _mm_slli_epi32(rcon, 1);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0xAB, 0xAB, 0xAB, 0xAB);
    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0x4D, 0x4D, 0x4D, 0x4D);
    {
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        ((__m128i *)output)[0] = _mm_aesenclast_si128(tmp0_0, round_key[0]);
        ((__m128i *)output)[1] = _mm_aesenclast_si128(tmp0_1, round_key[1]);
        ((__m128i *)output)[2] = _mm_aesenclast_si128(tmp1_0, round_key[0]);
        ((__m128i *)output)[3] = _mm_aesenclast_si128(tmp1_1, round_key[1]);
    }
}

static inline void rijndael_commit(uint8_t dst[2][32], const uint8_t salt[32], const uint32_t idx, const uint8_t seed[32]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    __m128i block_0[2] = {0};
    __m128i block_1[2] = {0};

    block256_t *const output = dst;
    const uint8_t *const seed_with_zeros = seed;

    uint8_t *msg = (uint8_t *)block_0;

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 32);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;

    // salt ^ (domain_separator || idx || 1)
    block_1[0] = block_0[0];
    block_1[1] = block_0[1];
    msg = (uint8_t *)block_1;
    msg[0] ^= 0x01;

    __m128i round_key[2];
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);
    __m128i RIJNDAEL256_MASK = _mm_set_epi32(0x03020d0c, 0x0f0e0908, 0x0b0a0504, 0x07060100);
    __m128i BLEND_MASK = _mm_set_epi32(0x80000000, 0x80800000, 0x80800000, 0x80808000);

    round_key[0] = ((__m128i *)seed_with_zeros)[0];
    round_key[1] = ((__m128i *)seed_with_zeros)[1];

    block_0[0] = _mm_xor_si128(block_0[0], round_key[0]); /* round 0 (initial xor) */
    block_0[1] = _mm_xor_si128(block_0[1], round_key[1]);
    block_1[0] = _mm_xor_si128(block_1[0], round_key[0]); /* round 0 (initial xor) */
    block_1[1] = _mm_xor_si128(block_1[1], round_key[1]);

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    for (int i = 1; i < 9; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);
        rcon = _mm_slli_epi32(rcon, 1);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0x1B, 0x1B, 0x1B, 0x1B);
    for (int i = 9; i < 13; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);
        rcon = _mm_slli_epi32(rcon, 1);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0xAB, 0xAB, 0xAB, 0xAB);
    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        block_0[0] = _mm_aesenc_si128(tmp0_0, round_key[0]);
        block_0[1] = _mm_aesenc_si128(tmp0_1, round_key[1]);
        block_1[0] = _mm_aesenc_si128(tmp1_0, round_key[0]);
        block_1[1] = _mm_aesenc_si128(tmp1_1, round_key[1]);
    }

    rcon = _mm_set_epi32(0x4D, 0x4D, 0x4D, 0x4D);
    {
        __m128i tmp = _mm_shuffle_epi8(round_key[1], shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        rijndael_256_assist(round_key, tmp, round_key);

        __m128i tmp0_0 = _mm_blendv_epi8(block_0[0], block_0[1], BLEND_MASK);
        __m128i tmp0_1 = _mm_blendv_epi8(block_0[1], block_0[0], BLEND_MASK);
        tmp0_0 = _mm_shuffle_epi8(tmp0_0, RIJNDAEL256_MASK);
        tmp0_1 = _mm_shuffle_epi8(tmp0_1, RIJNDAEL256_MASK);

        __m128i tmp1_0 = _mm_blendv_epi8(block_1[0], block_1[1], BLEND_MASK);
        __m128i tmp1_1 = _mm_blendv_epi8(block_1[1], block_1[0], BLEND_MASK);
        tmp1_0 = _mm_shuffle_epi8(tmp1_0, RIJNDAEL256_MASK);
        tmp1_1 = _mm_shuffle_epi8(tmp1_1, RIJNDAEL256_MASK);

        ((__m128i *)output)[0] = _mm_aesenclast_si128(tmp0_0, round_key[0]);
        ((__m128i *)output)[1] = _mm_aesenclast_si128(tmp0_1, round_key[1]);
        ((__m128i *)output)[2] = _mm_aesenclast_si128(tmp1_0, round_key[0]);
        ((__m128i *)output)[3] = _mm_aesenclast_si128(tmp1_1, round_key[1]);
    }
}

static inline void rijndael_expand_share(uint8_t (*dst)[32], const uint8_t salt[32], const uint8_t seed[32], uint8_t len) {
    // This function assumes dst has capacity len, and that the len is at most 255

    rijndael_256_round_keys_t key = {0};
    block256_t ctr;
    block256_t seed_vec;
    block256_t salt_vec;

    *((__m256i*)ctr) = _mm256_setzero_si256();
    memcpy(seed_vec, seed, sizeof(block256_t));
    memcpy(salt_vec, salt, sizeof(block256_t));

    rijndael_256_key_expansion(&key, seed_vec);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;
        block256_t msg;
        *((__m256i*)msg) = _mm256_setzero_si256();
        ((__m256i *)msg)[0] = _mm256_xor_si256(*(__m256i *)ctr, ((__m256i *)salt_vec)[0]);

        block256_t output;
        *((__m256i*)output) = _mm256_setzero_si256();
        rijndael_256_encrypt(output, msg, &key);
        memcpy(dst[i], output, sizeof(uint8_t) * 32);
    }
}


#endif //MIRATH_SEED_EXPAND_RIJNDAEL_256_H
