#ifndef HASH_SHA3_H
#define HASH_SHA3_H

#include <stdint.h>
#include "KeccakHash.h"

typedef Keccak_HashInstance hash_sha3_ctx;
typedef Keccak_HashInstance shake_prng_t;

#include "KeccakHashtimes4.h"

typedef Keccak_HashInstancetimes4 shake_prng_x4_t;
typedef Keccak_HashInstancetimes4 shake_sha3_x4_t;

void hash_sha3_x4_init(shake_sha3_x4_t *ctx);
void hash_sha3_x4_update(shake_sha3_x4_t *ctx, const uint8_t **input, size_t size);
void hash_sha3_x4_finalize(uint8_t **output, shake_sha3_x4_t *ctx);

void seedexpander_shake_x4_init(shake_prng_x4_t* seedexpander_shake, const uint8_t** seed, size_t seed_size, const uint8_t** salt, size_t salt_size);
void seedexpander_shake_x4_get_bytes(shake_prng_x4_t* seedexpander_shake, uint8_t** output, size_t output_size);

void seedexpander_shake_init(shake_prng_t* seedexpander_shake, const uint8_t* seed, size_t seed_size, const uint8_t* salt, size_t salt_size);
void seedexpander_shake_get_bytes(shake_prng_t* seedexpander_shake, uint8_t* output, size_t output_size);

void shake_init(shake_prng_t *prng);
void hash_init(hash_sha3_ctx *ctx);
void hash_update(hash_sha3_ctx *ctx, const uint8_t *input, uint32_t size);
void hash_finalize(uint8_t *output, hash_sha3_ctx *ctx);
void hash_shake(uint8_t *output, uint32_t output_size, const uint8_t * input, uint32_t input_size);
void hash_squeeze(hash_sha3_ctx *prng, void *target, uint32_t length);

#endif // HASH_SHA3_H
