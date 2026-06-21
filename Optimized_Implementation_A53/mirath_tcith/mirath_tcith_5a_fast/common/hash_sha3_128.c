/**
 * @file hash_sha3_128.c
 * @brief SHA3/SHAKE-128 hash functions for security level 1 parameter sets.
 *
 * Includes x4 sequential wrappers for API compatibility with AVX2 optimized version.
 * On ARM, x4 operations execute as 4 sequential single-instance calls.
 */

#include <string.h>
#include "hash_sha3.h"

/* ============================================================
 * Single-instance API (unchanged from reference implementation)
 * ============================================================ */

void seedexpander_shake_init(shake_prng_t* seedexpander_shake, const uint8_t* seed, size_t seed_size, const uint8_t* salt, size_t salt_size) {
    Keccak_HashInitialize_SHAKE128(seedexpander_shake);
    Keccak_HashUpdate(seedexpander_shake, seed, seed_size << 3);
    if (salt != NULL) Keccak_HashUpdate(seedexpander_shake, salt, salt_size << 3);
    Keccak_HashFinal(seedexpander_shake, NULL);
}

void seedexpander_shake_get_bytes(shake_prng_t* seedexpander_shake, uint8_t* output, size_t output_size) {
    Keccak_HashSqueeze(seedexpander_shake, output, output_size << 3);
}

void hash_init(hash_sha3_ctx *ctx) {
    Keccak_HashInitialize_SHA3_256(ctx);
}

void hash_update(hash_sha3_ctx *ctx, const uint8_t *input, uint32_t size) {
    Keccak_HashUpdate(ctx, input, size << 3);
}

void hash_finalize(uint8_t *output, hash_sha3_ctx *ctx) {
    Keccak_HashFinal(ctx, output);
}

void hash_shake(uint8_t *output, uint32_t output_size, const uint8_t *input, uint32_t input_size) {
    shake_prng_t ctx;
    Keccak_HashInitialize_SHAKE128(&ctx);
    Keccak_HashUpdate(&ctx, input, input_size << 3);
    Keccak_HashFinal(&ctx, NULL);
    Keccak_HashSqueeze(&ctx, output, output_size << 3);
}

void hash_squeeze(hash_sha3_ctx *prng, void *target, uint32_t length) {
    Keccak_HashSqueeze(prng, target, length << 3);
}

/* ============================================================
 * x4 Sequential Wrappers (ARM replacement for AVX2 parallel)
 *
 * These provide the same API as KeccakHashtimes4 but execute
 * as 4 sequential single-instance calls. This preserves the
 * optimized TCitH core code structure.
 * ============================================================ */

void hash_sha3_x4_init(shake_sha3_x4_t *ctx) {
    for (int i = 0; i < 4; i++) {
        Keccak_HashInitialize_SHA3_256(&ctx->instances[i]);
    }
}

void hash_sha3_x4_update(shake_sha3_x4_t *ctx, const uint8_t **input, size_t size) {
    for (int i = 0; i < 4; i++) {
        Keccak_HashUpdate(&ctx->instances[i], input[i], size << 3);
    }
}

void hash_sha3_x4_finalize(uint8_t **output, shake_sha3_x4_t *ctx) {
    for (int i = 0; i < 4; i++) {
        Keccak_HashFinal(&ctx->instances[i], output ? output[i] : NULL);
    }
}

void seedexpander_shake_x4_init(shake_prng_x4_t* seedexpander_shake, const uint8_t** seed, size_t seed_size, const uint8_t** salt, size_t salt_size) {
    for (int i = 0; i < 4; i++) {
        Keccak_HashInitialize_SHAKE128(&seedexpander_shake->instances[i]);
        Keccak_HashUpdate(&seedexpander_shake->instances[i], seed[i], seed_size << 3);
        if (salt != NULL) {
            Keccak_HashUpdate(&seedexpander_shake->instances[i], salt[i], salt_size << 3);
        }
        Keccak_HashFinal(&seedexpander_shake->instances[i], NULL);
    }
}

void seedexpander_shake_x4_get_bytes(shake_prng_x4_t* seedexpander_shake, uint8_t** output, size_t output_size) {
    for (int i = 0; i < 4; i++) {
        Keccak_HashSqueeze(&seedexpander_shake->instances[i], output[i], output_size << 3);
    }
}
