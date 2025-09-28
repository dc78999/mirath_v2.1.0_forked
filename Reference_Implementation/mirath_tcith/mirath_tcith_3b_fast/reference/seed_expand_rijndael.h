
#ifndef MIRATH_SEED_EXPAND_RIJNDAEL_192_H
#define MIRATH_SEED_EXPAND_RIJNDAEL_192_H

#include "rijndael.h"
#define DOMAIN_SEPARATOR_PRG 4
#define DOMAIN_SEPARATOR_CMT 3

static inline void rijndael_expand_seed(uint8_t dst[2][24], const uint8_t salt[24], const uint32_t idx, const uint8_t seed[24]) {
    // We assume that the output dst contains zeros

    uint8_t output[2][32] = {0};

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_PRG;
    uint8_t msg[32] = {0};
    aes_round_keys_t key;

    uint8_t seed_with_zeros[32] = {0};
    memcpy(seed_with_zeros, seed, sizeof(uint8_t) * 24);    // key = (seed || 0)

    rijndael_256_key_expansion(&key, seed_with_zeros);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 24);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    rijndael_256_encrypt(&key, msg, output[0]);
    memcpy(dst[0], output[0], sizeof(uint8_t) * 24); // copy only 24 bytes of output[0]

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    rijndael_256_encrypt(&key, msg, output[1]);
    memcpy(dst[1], output[1], sizeof(uint8_t) * 24); // copy only 24 bytes of output[1]

}

static inline void rijndael_commit(uint8_t dst[2][24], const uint8_t salt[24], const uint32_t idx, const uint8_t seed[24]) {
    // We assume that the output dst contains zeros

    uint8_t output[2][32] = {0};

    uint8_t domain_separator = (uint8_t)DOMAIN_SEPARATOR_CMT;
    uint8_t msg[32] = {0};
    aes_round_keys_t key;

    uint8_t seed_with_zeros[32] = {0};
    memcpy(seed_with_zeros, seed, sizeof(uint8_t) * 24);    // key = (seed || 0)

    rijndael_256_key_expansion(&key, seed_with_zeros);

    // salt ^ (domain_separator || idx || 0)
    memcpy(msg, salt, sizeof(uint8_t) * 24);
    msg[0] ^= 0x00;
    for (size_t k = 0; k < 4; k++) {
        msg[k + 1] ^= ((uint8_t *)&idx)[k];
    }
    msg[5] ^= domain_separator;
    rijndael_256_encrypt(&key, msg, output[0]);
    memcpy(dst[0], output[0], sizeof(uint8_t) * 24); // copy only 24 bytes of output[0]

    // salt ^ (domain_separator || idx || 1)
    msg[0] ^= 0x01;
    rijndael_256_encrypt(&key, msg, output[1]);
    memcpy(dst[1], output[1], sizeof(uint8_t) * 24); // copy only 24 bytes of output[1]

}

static inline void rijndael_expand_share(uint8_t (*dst)[24], const uint8_t salt[24], const uint8_t seed[24], uint8_t len) {
    // This function assumes dst has capacity len, and that the len is at most 255

    aes_round_keys_t key;
    uint8_t ctr[32] = {0};

    uint8_t seed_with_zeros[32] = {0};
    memcpy(seed_with_zeros, seed, sizeof(uint8_t) * 24);    // key = (0 || seed)

    rijndael_256_key_expansion(&key, seed_with_zeros);

    for (uint8_t i = 0; i < len; i++) {
        ctr[0] = i;

        uint8_t msg[32] = {0};
        for (size_t k = 0; k < 24; k++) {
            msg[k] = (ctr[k] ^ salt[k]);
        }

        uint8_t output[32] = {0};
        rijndael_256_encrypt(&key, msg, output);
        memcpy(dst[i], output, sizeof(uint8_t) * 24); // copy only 24 bytes of output
    }
}

#endif //MIRATH_SEED_EXPAND_RIJNDAEL_192_H
