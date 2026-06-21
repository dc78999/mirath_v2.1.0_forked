#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "api.h"
#include "bench_cycles.h"
#include "rng.h"

#define N_ITER 100

static int cmp_uint64(const void *a, const void *b) {
    if (*(uint64_t*)a < *(uint64_t*)b) return -1;
    if (*(uint64_t*)a > *(uint64_t*)b) return 1;
    return 0;
}

int main(void) {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char m[32] = {0};
    unsigned char sm[CRYPTO_BYTES + 32];
    unsigned long long smlen, mlen;
    unsigned char seed[48] = {0};
    
    uint64_t t1, t2;
    uint64_t cycles_keygen[N_ITER];
    uint64_t cycles_sign[N_ITER];
    uint64_t cycles_verify[N_ITER];

    randombytes_init(seed, NULL, 256);

    printf("Benchmarking Mirath (N=%d)\n", N_ITER);

    for (int i = 0; i < N_ITER; i++) {
        t1 = cpucyclesStart();
        if (crypto_sign_keypair(pk, sk) != 0) {
            printf("Keygen failed\n");
            return -1;
        }
        t2 = cpucyclesStop();
        cycles_keygen[i] = t2 - t1;

        t1 = cpucyclesStart();
        if (crypto_sign(sm, &smlen, m, 32, sk) != 0) {
            printf("Sign failed\n");
            return -1;
        }
        t2 = cpucyclesStop();
        cycles_sign[i] = t2 - t1;

        t1 = cpucyclesStart();
        if (crypto_sign_open(m, &mlen, sm, smlen, pk) != 0) {
            printf("Verify failed\n");
            return -1;
        }
        t2 = cpucyclesStop();
        cycles_verify[i] = t2 - t1;
    }

    qsort(cycles_keygen, N_ITER, sizeof(uint64_t), cmp_uint64);
    qsort(cycles_sign, N_ITER, sizeof(uint64_t), cmp_uint64);
    qsort(cycles_verify, N_ITER, sizeof(uint64_t), cmp_uint64);

    printf("Median KeyGen: %" PRIu64 " cycles\n", cycles_keygen[N_ITER/2]);
    printf("Median Sign:   %" PRIu64 " cycles\n", cycles_sign[N_ITER/2]);
    printf("Median Verify: %" PRIu64 " cycles\n", cycles_verify[N_ITER/2]);

    return 0;
}
