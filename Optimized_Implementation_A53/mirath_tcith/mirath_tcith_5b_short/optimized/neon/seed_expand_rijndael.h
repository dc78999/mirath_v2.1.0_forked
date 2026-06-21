/**
 * @file seed_expand_rijndael.h
 * @brief Seed expand functions based on AES-128 using ARM NEON Cryptography extensions
 */

#ifndef SEED_EXPAND_RIJNDAEL_128_H
#define SEED_EXPAND_RIJNDAEL_128_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

typedef uint8_t block128_t[16] __attribute__ ((aligned (16)));

static inline void rijndael_expand_seed(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    uint8x16_t block_0 = vdupq_n_u8(0);
    uint8x16_t block_1 = vdupq_n_u8(0);

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

    uint8x16_t round_key = vld1q_u8(seed);
    uint8x16_t prev_key = round_key;

    const uint8_t shuffle_arr[16] = {
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c
    };
    uint8x16_t shuffle_mask = vld1q_u8(shuffle_arr);
    uint8x16_t rcon = vreinterpretq_u8_u32(vdupq_n_u32(1));

    for (int i = 1; i < 9; i++) {
        block_0 = vaesmcq_u8(vaeseq_u8(block_0, prev_key));
        block_1 = vaesmcq_u8(vaeseq_u8(block_1, prev_key));

        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
        
        prev_key = round_key;
    }

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0x1B));
    {
        block_0 = vaesmcq_u8(vaeseq_u8(block_0, prev_key));
        block_1 = vaesmcq_u8(vaeseq_u8(block_1, prev_key));

        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
        
        prev_key = round_key;
    }

    {
        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);

        uint8x16_t out0 = veorq_u8(vaeseq_u8(block_0, prev_key), round_key);
        uint8x16_t out1 = veorq_u8(vaeseq_u8(block_1, prev_key), round_key);

        vst1q_u8(dst[0], out0);
        vst1q_u8(dst[1], out1);
    }
}

static inline void rijndael_commit(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    uint8x16_t block_0 = vdupq_n_u8(0);
    uint8x16_t block_1 = vdupq_n_u8(0);

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

    uint8x16_t round_key = vld1q_u8(seed);
    uint8x16_t prev_key = round_key;

    const uint8_t shuffle_arr[16] = {
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c
    };
    uint8x16_t shuffle_mask = vld1q_u8(shuffle_arr);
    uint8x16_t rcon = vreinterpretq_u8_u32(vdupq_n_u32(1));

    for (int i = 1; i < 9; i++) {
        block_0 = vaesmcq_u8(vaeseq_u8(block_0, prev_key));
        block_1 = vaesmcq_u8(vaeseq_u8(block_1, prev_key));

        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
        
        prev_key = round_key;
    }

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0x1B));
    {
        block_0 = vaesmcq_u8(vaeseq_u8(block_0, prev_key));
        block_1 = vaesmcq_u8(vaeseq_u8(block_1, prev_key));

        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);
        rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
        
        prev_key = round_key;
    }

    {
        uint8x16_t tmp = vqtbl1q_u8(round_key, shuffle_mask);
        tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
        round_key = aes_128_assist(round_key, tmp);

        uint8x16_t out0 = veorq_u8(vaeseq_u8(block_0, prev_key), round_key);
        uint8x16_t out1 = veorq_u8(vaeseq_u8(block_1, prev_key), round_key);

        vst1q_u8(dst[0], out0);
        vst1q_u8(dst[1], out1);
    }
}

static inline void rijndael_expand_share(uint8_t (*dst)[16], const uint8_t salt[16], const uint8_t seed[16], uint8_t len) {
    aes_128_round_keys_t key = {0};
    block128_t ctr;
    block128_t seed_vec;
    block128_t salt_vec;

    memset(ctr, 0, sizeof(block128_t));
    memcpy(seed_vec, seed, sizeof(block128_t));
    memcpy(salt_vec, salt, sizeof(block128_t));

    aes_128_key_expansion(&key, seed_vec);

    uint8x16_t v_salt = vld1q_u8(salt_vec);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;
        uint8x16_t v_ctr = vld1q_u8(ctr);
        uint8x16_t msg = veorq_u8(v_ctr, v_salt);

        block128_t output;
        aes_128_encrypt(output, (uint8_t *)&msg, &key);
        memcpy(dst[i], output, sizeof(uint8_t) * 16);
    }
}

#endif //SEED_EXPAND_RIJNDAEL_128_H
