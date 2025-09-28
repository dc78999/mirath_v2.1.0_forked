/**
 * @file seed_expand_functions_ref.h
 * @brief Content for seed_expand_functions_ref.h (Seed expand functions based on AES-128 and Rijndael-256)
 */

#ifndef MIRATH_SEED_EXPAND_RIJNDAEL_128_H
#define MIRATH_SEED_EXPAND_RIJNDAEL_128_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

static inline void rijndael_expand_seed(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    uint8_t msg[16] = {0};
    aes_round_keys_t key;

    aes_128_key_expansion(&key, seed);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 16);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    aes_128_encrypt(&key, msg, dst[0]);

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    aes_128_encrypt(&key, msg, dst[1]);

}

static inline void rijndael_commit(uint8_t dst[2][16], const uint8_t salt[16], const uint32_t idx, const uint8_t seed[16]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    uint8_t msg[16] = {0};
    aes_round_keys_t key;

    aes_128_key_expansion(&key, seed);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 16);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    aes_128_encrypt(&key, msg, dst[0]);

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    aes_128_encrypt(&key, msg, dst[1]);

}

static inline void rijndael_expand_share(uint8_t (*dst)[16], const uint8_t salt[16], const uint8_t seed[16], uint8_t len) {
    // This function assumes dst has capacity len, and that the len is at most 255

    aes_round_keys_t key;
    uint8_t ctr[16] = {0};

    aes_128_key_expansion(&key, seed);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;

        uint8_t msg[16] = {0};
        for (size_t k = 0; k < 16; k++) {
            msg[k] = (ctr[k] ^ salt[k]);
        }

        aes_128_encrypt(&key, msg, dst[i]);
    }
}

#endif //MIRATH_SEED_EXPAND_RIJNDAEL_128_H
