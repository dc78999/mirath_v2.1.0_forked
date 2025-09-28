/**
 * @file seed_expand_functions_avx.h
 * @brief Content for seed_expand_functions_avx.h (Seed expand functions based on AES-128)
 */

#ifndef SEED_EXPAND_RIJNDAEL_128_H
#define SEED_EXPAND_RIJNDAEL_128_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

typedef uint8_t block128_t[16] __attribute__ ((aligned (16)));

static inline void rijndael_expand_seed(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    __m128i block_0 = {0};
    __m128i block_1 = {0};

    uint8_t *msg = (uint8_t *)&block_0;

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 16);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;

    // salt ^ (domain_separator || idx || 1)
    block_1 = block_0;
    msg = (uint8_t *)&block_1;
    msg[0] ^= 0x01;

    __m128i round_key;
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);

    round_key = *((__m128i *)seed);
    block_0 = _mm_xor_si128(block_0, round_key); /* round 0 (initial xor) */
    block_1 = _mm_xor_si128(block_1, round_key); /* round 0 (initial xor) */

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    for (int i = 1; i < 9; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = _mm_slli_epi32(rcon, 1);

        block_0 = _mm_aesenc_si128(block_0, round_key);
        block_1 = _mm_aesenc_si128(block_1, round_key);
    }

    rcon = _mm_set_epi32(0x1B, 0x1B, 0x1B, 0x1B);
    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = _mm_slli_epi32(rcon, 1);

        block_0 = _mm_aesenc_si128(block_0, round_key);
        block_1 = _mm_aesenc_si128(block_1, round_key);
    }

    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);

        ((__m128i *)dst)[0] = _mm_aesenclast_si128(block_0, round_key);
        ((__m128i *)dst)[1] = _mm_aesenclast_si128(block_1, round_key);
    }
}

static inline void rijndael_commit(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    __m128i block_0 = {0};
    __m128i block_1 = {0};

    uint8_t *msg = (uint8_t *)&block_0;

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 16);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;

    // salt ^ (domain_separator || idx || 1)
    block_1 = block_0;
    msg = (uint8_t *)&block_1;
    msg[0] ^= 0x01;

    __m128i round_key;
    __m128i shuffle_mask = _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d);

    round_key = *((__m128i *)seed);
    block_0 = _mm_xor_si128(block_0, round_key); /* round 0 (initial xor) */
    block_1 = _mm_xor_si128(block_1, round_key); /* round 0 (initial xor) */

    __m128i rcon = _mm_set_epi32(1, 1, 1, 1);
    for (int i = 1; i < 9; i++) {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = _mm_slli_epi32(rcon, 1);

        block_0 = _mm_aesenc_si128(block_0, round_key);
        block_1 = _mm_aesenc_si128(block_1, round_key);
    }

    rcon = _mm_set_epi32(0x1B, 0x1B, 0x1B, 0x1B);
    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = _mm_slli_epi32(rcon, 1);

        block_0 = _mm_aesenc_si128(block_0, round_key);
        block_1 = _mm_aesenc_si128(block_1, round_key);
    }

    {
        // on the fly key scheduling
        __m128i tmp = _mm_shuffle_epi8(round_key, shuffle_mask);
        tmp = _mm_aesenclast_si128(tmp, rcon);
        round_key = aes_128_assist(round_key, tmp);

        ((__m128i *)dst)[0] = _mm_aesenclast_si128(block_0, round_key);
        ((__m128i *)dst)[1] = _mm_aesenclast_si128(block_1, round_key);
    }
}

static inline void rijndael_expand_share(uint8_t (*dst)[16], const uint8_t salt[16], const uint8_t seed[16], uint8_t len) {
    // This function assumes dst has capacity len, and that the len is at most 255

    aes_128_round_keys_t key = {0};
    block128_t ctr;
    block128_t seed_vec;
    block128_t salt_vec;

    *((__m128i*)ctr) = _mm_setzero_si128();
    memcpy(seed_vec, seed, sizeof(block128_t));
    memcpy(salt_vec, salt, sizeof(block128_t));

    aes_128_key_expansion(&key, seed_vec);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;
        block128_t msg;
        *((__m128i*)msg) = _mm_setzero_si128();
        ((__m128i *)msg)[0] = _mm_xor_si128(*(__m128i *)ctr, ((__m128i *)salt_vec)[0]);

        block128_t output;
        *((__m128i*)output) = _mm_setzero_si128();
        aes_128_encrypt(output, msg, &key);
        memcpy(dst[i], output, sizeof(uint8_t) * 16);
    }
}

#endif //SEED_EXPAND_RIJNDAEL_128_H
