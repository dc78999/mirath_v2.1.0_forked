#include "hash_sha3.h"

void shake_init(shake_prng_t *prng) {
    Keccak_HashInitialize_SHAKE256(prng);
}

void hash_init(hash_sha3_ctx *ctx) {
    Keccak_HashInitialize_SHA3_512(ctx);
}

void hash_update(hash_sha3_ctx *ctx, const uint8_t *input, uint32_t size) {
    Keccak_HashUpdate(ctx, input, size << 3);
}

void hash_finalize(uint8_t *output, hash_sha3_ctx *ctx) {
    Keccak_HashFinal(ctx, output);
}

void hash_shake(uint8_t *output, uint32_t output_size, const uint8_t * input, uint32_t input_size) {
    SHAKE256(output, output_size, input, input_size);
}

void seedexpander_shake_init(shake_prng_t* seedexpander_shake, const uint8_t* seed, size_t seed_size, const uint8_t* salt, size_t salt_size) {
    Keccak_HashInitialize_SHAKE256(seedexpander_shake);
    Keccak_HashUpdate(seedexpander_shake, seed, seed_size << 3);
    if(salt != NULL) Keccak_HashUpdate(seedexpander_shake, salt, salt_size << 3);
    Keccak_HashFinal(seedexpander_shake, NULL);
}

void seedexpander_shake_get_bytes(shake_prng_t* seedexpander_shake, uint8_t* output, size_t output_size) {
    Keccak_HashSqueeze(seedexpander_shake, output, output_size << 3);
}

