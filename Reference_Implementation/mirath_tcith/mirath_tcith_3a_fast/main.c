#define _GNU_SOURCE

#include <unistd.h>
#include <sys/syscall.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "rng.h"
#include "api.h"

inline static uint64_t cpucyclesStart (void) {
    unsigned hi, lo;
    __asm__ __volatile__ (	"CPUID\n\t"
                              "RDTSC\n\t"
                              "mov %%edx, %0\n\t"
                              "mov %%eax, %1\n\t"
            : "=r" (hi), "=r" (lo)
            :
            : "%rax", "%rbx", "%rcx", "%rdx");

    return ((uint64_t) lo) ^ (((uint64_t) hi) << 32);
}



inline static uint64_t cpucyclesStop (void) {
    unsigned hi, lo;
    __asm__ __volatile__(	"RDTSCP\n\t"
                             "mov %%edx, %0\n\t"
                             "mov %%eax, %1\n\t"
                             "CPUID\n\t"
            : "=r" (hi), "=r" (lo)
            :
            : "%rax", "%rbx", "%rcx", "%rdx");

    return ((uint64_t) lo) ^ (((uint64_t) hi) << 32);
}

int main(void) {

    unsigned long long mlen = 22;
    unsigned long long smlen;

    unsigned char m[] = {0x52, 0x61, 0x6e, 0x6b, 0x20, 0x73, 0x59, 0x6e, 0x64, 0x72, 0x6f, 0x6d,
                        0x65, 0x20, 0x44, 0x45, 0x63, 0x6f, 0x64, 0x69, 0x6e, 0x67};

    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char sm[CRYPTO_BYTES + mlen];

    unsigned long long t1, t2, t3, t4, t5, t6;

    unsigned char seed[48] = {0};
//    (void)syscall(SYS_getrandom, seed, 48, 0);
    randombytes_init(seed, NULL, 256);

    /***************/
    /*   Mirath    */
    /***************/

    t1 = cpucyclesStart();
    if (crypto_sign_keypair(pk, sk) != 0) {
        printf("crypto_sign_keypair: Failed\n");
        return -1;
    }
    t2 = cpucyclesStop();

    t3 = cpucyclesStart();
    if (crypto_sign(sm, &smlen, m, mlen, sk) != 0) {
        printf("crypto_sign: Failed\n");
        return -1;
    }
    t4 = cpucyclesStart();

    t5 = cpucyclesStart();
    if (crypto_sign_open(m, &mlen, sm, smlen, pk) != 0) {
        printf("crypto_sign_open: Failed\n");
        return -1;
    }
    t6 = cpucyclesStart();

    printf("\n mirath_tcith_3a_fast ");
    printf("\n  crypto_sign_keypair: %" PRIu64 " CPU cycles", (uint64_t)(t2 - t1));
    printf("\n  crypto_sign:         %" PRIu64 " CPU cycles", (uint64_t)(t4 - t3));
    printf("\n  crypto_sign_open:    %" PRIu64 " CPU cycles", (uint64_t)(t6 - t5));
    printf("\n\n");
    printf("\n sk: "); for(int k = 0 ; k < CRYPTO_SECRETKEYBYTES ; ++k) printf("%02x", sk[k]);
    printf("\n pk: "); for(int k = 0 ; k < CRYPTO_PUBLICKEYBYTES ; ++k) printf("%02x", pk[k]);
    printf("\n  m: "); for(int k = 0 ; k < (int)mlen ; ++k) printf("%02x", m[k]);
    printf("\n sm: "); for(int k = 0 ; k < ((int)smlen) ; ++k) printf("%02x", sm[k]);

    printf("\n\n");
    return 0;
}

