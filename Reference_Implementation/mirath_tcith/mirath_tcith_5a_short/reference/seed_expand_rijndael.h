
#ifndef MIRATH_SEED_EXPAND_RIJNDAEL_256_H
#define MIRATH_SEED_EXPAND_RIJNDAEL_256_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

static inline void rijndael_expand_seed(uint8_t dst[2][32], const uint8_t salt[32], const uint32_t idx, const uint8_t seed[32]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    uint8_t msg[32] = {0};
    aes_round_keys_t key;

    rijndael_256_key_expansion(&key, seed);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 32);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    rijndael_256_encrypt(&key, msg, dst[0]);

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    rijndael_256_encrypt(&key, msg, dst[1]);

}

static inline void rijndael_commit(uint8_t dst[2][32], const uint8_t salt[32], const uint32_t idx, const uint8_t seed[32]) {
    // We assume that the output dst contains zeros

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    uint8_t msg[32] = {0};
    aes_round_keys_t key;

    rijndael_256_key_expansion(&key, seed);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 32);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    rijndael_256_encrypt(&key, msg, dst[0]);

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    rijndael_256_encrypt(&key, msg, dst[1]);

}

static inline void rijndael_expand_share(uint8_t (*dst)[32], const uint8_t salt[32], const uint8_t seed[32], uint8_t len) {
    // This function assumes dst has capacity len, and that the len is at most 255

    aes_round_keys_t key;
    uint8_t ctr[32] = {0};

    rijndael_256_key_expansion(&key, seed);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;

        uint8_t msg[32] = {0};
        for (size_t k = 0; k < 32; k++) {
            msg[k] = (ctr[k] ^ salt[k]);
        }

        rijndael_256_encrypt(&key, msg, dst[i]);
    }
}

#endif //MIRATH_SEED_EXPAND_RIJNDAEL_256_H
