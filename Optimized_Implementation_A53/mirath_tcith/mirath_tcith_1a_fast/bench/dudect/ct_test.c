#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"
#include "dudect.h"
#include "../bench_cycles.h"
#include "rng.h"

#define MEASUREMENTS 10000

int main(void) {
    printf("Running Dudect Side-Channel Evaluation (%d measurements)...\n", MEASUREMENTS);

    dudect_ctx_t ctx;
    dudect_init(&ctx);

    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char m[32] = {0};
    unsigned char sm[CRYPTO_BYTES + 32];
    unsigned long long smlen;
    unsigned char seed[48] = {0};

    randombytes_init(seed, NULL, 256);
    crypto_sign_keypair(pk, sk);

    // Two classes of inputs: 
    // Class 0: random messages
    // Class 1: fixed message
    
    unsigned char class_0_msg[32];
    unsigned char class_1_msg[32];
    memset(class_1_msg, 0x42, 32);

    for (int i = 0; i < MEASUREMENTS; i++) {
        uint8_t class_label = i % 2;
        randombytes(class_0_msg, 32);
        
        unsigned char *test_msg = class_label ? class_1_msg : class_0_msg;

        uint64_t t1 = cpucyclesStart();
        crypto_sign(sm, &smlen, test_msg, 32, sk);
        uint64_t t2 = cpucyclesStop();

        dudect_push(&ctx, class_label, (double)(t2 - t1));
    }

    int ret = dudect_report(&ctx);
    if (ret == 0) {
        printf("Constant-time check PASSED for crypto_sign.\n");
    } else {
        printf("Constant-time check FAILED for crypto_sign.\n");
    }

    return ret;
}
