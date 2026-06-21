#include <string.h>
#include "rng.h"
#include <stdio.h>
#include <stdlib.h>

int seedexpander_init(AES_XOF_struct *ctx,
                  unsigned char *seed,
                  unsigned char *diversifier,
                  unsigned long maxlen)
{
    (void)ctx; (void)seed; (void)diversifier; (void)maxlen;
    return RNG_SUCCESS;
}

int seedexpander(AES_XOF_struct *ctx, unsigned char *x, unsigned long xlen)
{
    (void)ctx; (void)x; (void)xlen;
    return RNG_SUCCESS;
}

void randombytes_init(unsigned char *entropy_input,
                 unsigned char *personalization_string,
                 int security_strength)
{
    (void)entropy_input; (void)personalization_string; (void)security_strength;
}

int randombytes(unsigned char *x, unsigned long long xlen)
{
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) {
        perror("fopen");
        exit(1);
    }
    if (fread(x, 1, xlen, f) != xlen) {
        fclose(f);
        return RNG_BAD_OUTBUF;
    }
    fclose(f);
    return RNG_SUCCESS;
}
