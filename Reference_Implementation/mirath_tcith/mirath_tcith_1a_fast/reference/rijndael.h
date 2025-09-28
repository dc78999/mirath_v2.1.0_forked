/**
 * @file rijndael_ref.h
 * @brief Content for rijndael_ref.h (AES-128 and Rijndael-256 reference implementation)
 */

#ifndef RIJNDAEL_REF_H
#define RIJNDAEL_REF_H

#include <stdint.h>
#include <stdlib.h>

#define AES_MAX_ROUNDS 14

typedef uint8_t aes_word_t[4];
typedef aes_word_t aes_round_key_t[8];

// # of rows
#define AES_NR 4

// block with 4 (AES) up to 8 (Rijndael-256) units
typedef aes_word_t aes_block_t[8];

typedef struct {
    aes_round_key_t round_keys[AES_MAX_ROUNDS + 1];
} aes_round_keys_t;

// begin fields

typedef uint8_t bf8_t;

// GF(2^8) with X^8 + X^4 + X^3 + X^1 + 1
#define bf8_modulus (UINT8_C((1 << 4) | (1 << 3) | (1 << 1) | 1))

static inline bf8_t bf8_load(const uint8_t* src) {
    return *src;
}

static inline void bf8_store(uint8_t* dst, bf8_t src) {
    *dst = src;
}

static inline bf8_t bf8_mul(bf8_t lhs, bf8_t rhs) {
    bf8_t result = -(rhs & 1) & lhs;
    for (unsigned int idx = 1; idx < 8; ++idx) {
        const uint8_t mask = -((lhs >> 7) & 1);
        lhs = (lhs << 1) ^ (mask & bf8_modulus);
        result ^= -((rhs >> idx) & 1) & lhs;
    }
    return result;
}

static inline bf8_t bf8_square(bf8_t lhs) {
    bf8_t result = -(lhs & 1) & lhs;
    bf8_t rhs = lhs;
    for (unsigned int idx = 1; idx < 8; ++idx) {
        const uint8_t mask = -((lhs >> 7) & 1);
        lhs                = (lhs << 1) ^ (mask & bf8_modulus);
        result ^= -((rhs >> idx) & 1) & lhs;
    }
    return result;
}

static inline bf8_t bf8_inv(bf8_t in) {
    const bf8_t t2 = bf8_square(in);
    const bf8_t t3 = bf8_mul(in, t2);
    const bf8_t t5 = bf8_mul(t3, t2);
    const bf8_t t7 = bf8_mul(t5, t2);
    const bf8_t t14 = bf8_square(t7);
    const bf8_t t28 = bf8_square(t14);
    const bf8_t t56 = bf8_square(t28);
    const bf8_t t63 = bf8_mul(t56, t7);
    const bf8_t t126 = bf8_square(t63);
    const bf8_t t252 = bf8_square(t126);
    return bf8_mul(t252, t2);
}

// end fields

static inline uint8_t parity8(uint8_t n) {
    n ^= n >> 4;
    n ^= n >> 2;
    n ^= n >> 1;
    return !((~n) & 1);
}

// rijndael

#define set_bit(value, index) ((value) << (index))

#include <string.h>

#define ROUNDS_128 10
#define ROUNDS_256 14

#define KEY_WORDS_128 4
#define KEY_WORDS_256 8

#define AES_BLOCK_WORDS 4
#define RIJNDAEL_BLOCK_WORDS_256 8

static inline void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        out[i] = a[i] ^ b[i];
    }
}

static const bf8_t round_constants[30] = {
        0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
        0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

static int contains_zero(const bf8_t* block) {
    return !block[0] | !block[1] | !block[2] | !block[3];
}

static bf8_t compute_sbox(bf8_t in) {
    bf8_t t = bf8_inv(in);
    bf8_t t0 = set_bit(parity8(t & (1 | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7))), 0);
    t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 5) | (1 << 6) | (1 << 7))), 1);
    t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 6) | (1 << 7))), 2);
    t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 7))), 3);
    t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))), 4);
    t0 ^= set_bit(parity8(t & ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5))), 5);
    t0 ^= set_bit(parity8(t & ((1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6))), 6);
    t0 ^= set_bit(parity8(t & ((1 << 3) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7))), 7);
    return t0 ^ (1 | (1 << 1) | (1 << 5) | (1 << 6));
}

// ## AES ##
// Round Functions
static void add_round_key(unsigned int round, aes_block_t state, const aes_round_keys_t* round_key,
                          unsigned int block_words) {
    for (unsigned int c = 0; c < block_words; c++) {
        xor_u8_array(&state[c][0], &round_key->round_keys[round][c][0], &state[c][0], AES_NR);
    }
}

static int sub_bytes(aes_block_t state, unsigned int block_words) {
    int ret = 0;

    for (unsigned int c = 0; c < block_words; c++) {
        ret |= contains_zero(&state[c][0]);
        for (unsigned int r = 0; r < AES_NR; r++) {
            state[c][r] = compute_sbox(state[c][r]);
        }
    }

    return ret;
}

static void shift_row(aes_block_t state, unsigned int block_words) {
    aes_block_t new_state;
    switch (block_words) {
        case 4:
        case 6:
            for (unsigned int i = 0; i < block_words; ++i) {
                new_state[i][0] = state[i][0];
                new_state[i][1] = state[(i + 1) % block_words][1];
                new_state[i][2] = state[(i + 2) % block_words][2];
                new_state[i][3] = state[(i + 3) % block_words][3];
            }
            break;
        case 8:
            for (unsigned int i = 0; i < block_words; i++) {
                new_state[i][0] = state[i][0];
                new_state[i][1] = state[(i + 1) % 8][1];
                new_state[i][2] = state[(i + 3) % 8][2];
                new_state[i][3] = state[(i + 4) % 8][3];
            }
            break;
    }

    for (unsigned int i = 0; i < block_words; ++i) {
        memcpy(&state[i][0], &new_state[i][0], AES_NR);
    }
}

static void mix_column(aes_block_t state, unsigned int block_words) {
    for (unsigned int c = 0; c < block_words; c++) {
        bf8_t tmp = bf8_mul(state[c][0], 0x02) ^ bf8_mul(state[c][1], 0x03) ^ state[c][2] ^ state[c][3];
        bf8_t tmp_1 =
                state[c][0] ^ bf8_mul(state[c][1], 0x02) ^ bf8_mul(state[c][2], 0x03) ^ state[c][3];
        bf8_t tmp_2 =
                state[c][0] ^ state[c][1] ^ bf8_mul(state[c][2], 0x02) ^ bf8_mul(state[c][3], 0x03);
        bf8_t tmp_3 =
                bf8_mul(state[c][0], 0x03) ^ state[c][1] ^ state[c][2] ^ bf8_mul(state[c][3], 0x02);

        state[c][0] = tmp;
        state[c][1] = tmp_1;
        state[c][2] = tmp_2;
        state[c][3] = tmp_3;
    }
}

// Key Expansion functions
static void sub_words(bf8_t* words) {
    words[0] = compute_sbox(words[0]);
    words[1] = compute_sbox(words[1]);
    words[2] = compute_sbox(words[2]);
    words[3] = compute_sbox(words[3]);
}

static void rot_word(bf8_t* words) {
    bf8_t tmp = words[0];
    words[0]  = words[1];
    words[1]  = words[2];
    words[2]  = words[3];
    words[3]  = tmp;
}

static inline int expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
               unsigned int block_words, unsigned int num_rounds) {
    int ret = 0;

    for (unsigned int k = 0; k < key_words; k++) {
        round_keys->round_keys[k / block_words][k % block_words][0] = bf8_load(&key[4 * k]);
        round_keys->round_keys[k / block_words][k % block_words][1] = bf8_load(&key[(4 * k) + 1]);
        round_keys->round_keys[k / block_words][k % block_words][2] = bf8_load(&key[(4 * k) + 2]);
        round_keys->round_keys[k / block_words][k % block_words][3] = bf8_load(&key[(4 * k) + 3]);
    }

    for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
        bf8_t tmp[AES_NR];
        memcpy(tmp, round_keys->round_keys[(k - 1) / block_words][(k - 1) % block_words], sizeof(tmp));

        if (k % key_words == 0) {
            rot_word(tmp);
            ret |= contains_zero(tmp);
            sub_words(tmp);
            tmp[0] ^= round_constants[(k / key_words) - 1];
        }

        if (key_words > 6 && (k % key_words) == 4) {
            ret |= contains_zero(tmp);
            sub_words(tmp);
        }

        unsigned int m = k - key_words;
        round_keys->round_keys[k / block_words][k % block_words][0] =
                round_keys->round_keys[m / block_words][m % block_words][0] ^ tmp[0];
        round_keys->round_keys[k / block_words][k % block_words][1] =
                round_keys->round_keys[m / block_words][m % block_words][1] ^ tmp[1];
        round_keys->round_keys[k / block_words][k % block_words][2] =
                round_keys->round_keys[m / block_words][m % block_words][2] ^ tmp[2];
        round_keys->round_keys[k / block_words][k % block_words][3] =
                round_keys->round_keys[m / block_words][m % block_words][3] ^ tmp[3];
    }

    return ret;
}

// Calling Functions

static inline int aes_128_key_expansion(aes_round_keys_t* round_key, const uint8_t* key) {
    return expand_key(round_key, key, KEY_WORDS_128, AES_BLOCK_WORDS, ROUNDS_128);
}

static inline int rijndael_256_key_expansion(aes_round_keys_t* round_key, const uint8_t* key) {
    return expand_key(round_key, key, KEY_WORDS_256, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
}

static void load_state(aes_block_t state, const uint8_t* src, unsigned int block_words) {
    for (unsigned int i = 0; i != block_words * 4; ++i) {
        state[i / 4][i % 4] = bf8_load(&src[i]);
    }
}

static void store_state(uint8_t* dst, aes_block_t state, unsigned int block_words) {
    for (unsigned int i = 0; i != block_words * 4; ++i) {
        bf8_store(&dst[i], state[i / 4][i % 4]);
    }
}

static int aes_encrypt(const aes_round_keys_t* keys, aes_block_t state, unsigned int block_words,
                       unsigned int num_rounds) {
    int ret = 0;

    // first round
    add_round_key(0, state, keys, block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
        ret |= sub_bytes(state, block_words);
        shift_row(state, block_words);
        mix_column(state, block_words);
        add_round_key(round, state, keys, block_words);
    }

    // last round
    ret |= sub_bytes(state, block_words);
    shift_row(state, block_words);
    add_round_key(num_rounds, state, keys, block_words);

    return ret;
}

static inline int aes_128_encrypt(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
    aes_block_t state;
    load_state(state, plaintext, AES_BLOCK_WORDS);
    const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_128);
    store_state(ciphertext, state, AES_BLOCK_WORDS);
    return ret;
}

static inline int rijndael_256_encrypt(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext) {
    aes_block_t state;
    load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_256);
    const int ret = aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
    store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_256);
    return ret;
}


#endif //RIJNDAEL_REF_H
