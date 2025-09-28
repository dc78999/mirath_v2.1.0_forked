#ifndef MIRATH_TCITH_API_H
#define MIRATH_TCITH_API_H

#include "mirath_parameters.h"

#define CRYPTO_ALGNAME "Mirath"

#define CRYPTO_PUBLICKEYBYTES MIRATH_PUBLIC_KEY_BYTES
#define CRYPTO_SECRETKEYBYTES MIRATH_SECRET_KEY_BYTES
#define CRYPTO_BYTES MIRATH_SIGNATURE_BYTES

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen, const unsigned char *sk);
int crypto_sign_open(unsigned char *msg, unsigned long long *mlen, const unsigned char *sig_msg, unsigned long long smlen, const unsigned char *pk);

#endif //MIRATH_TCITH_API_H
