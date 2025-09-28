/**
 * @file sign.c
 * @brief Implementation of the NIST api functions
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "api.h"
#include "mirath_keygen.h"
#include "mirath_sign.h"
#include "mirath_verify.h"

/**
 * @brief Generate a keypair.
 *
 * @param [out] pk pointer to public key bytes
 * @param [out] sk pointer to public key bytes
 * @returns 0 if key generation is successful
 */
int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {

    if (mirath_keygen(pk, sk) != EXIT_SUCCESS) {
        memset(sk, 0, CRYPTO_SECRETKEYBYTES);
        return -1;
    }
    return 0;
}

/**
 * @brief Generate a signature of a message.
 *
 * @param [out] sig_msg pointer to output signed message
 *                 (allocated array with CRYPTO_BYTES + mlen bytes), can be equal to m
 * @param [out] smlen pointer to output length of signed message
 * @param [in] msg pointer to message to be signed
 * @param [in] mlen length of the message
 * @param [in] sk pointer to the secret key bytes
 * @returns 0 if signing is successful
 */
int crypto_sign(unsigned char *sig_msg, unsigned long long *smlen, const unsigned char *msg, unsigned long long mlen, const unsigned char *sk) {
    if (mirath_sign(sig_msg, (uint8_t*)msg, mlen, (uint8_t*)sk) != EXIT_SUCCESS) {
        memset(sig_msg, 0, CRYPTO_BYTES + mlen);
        return -1;
    }

    for (size_t i = 0; i < mlen; ++i) { sig_msg[CRYPTO_BYTES + mlen - 1 - i] = msg[mlen - 1 - i]; }
    *smlen = mlen + CRYPTO_BYTES;
    return 0;
}

/**
 * @brief Verify a signed message
 *
 * @param [out] msg pointer to output message
 *                (allocated array with smlen bytes), can be equal to sm
 * @param [out] mlen pointer to output length of message
 * @param [in] sig_msg pointer to signed message
 * @param [in] smlen length of signed message
 * @param [in] pk pointer to the public key bytes
 * @returns 0 if signed message could be verified correctly and -1 otherwise
 */
int crypto_sign_open(unsigned char *msg, unsigned long long *mlen, const unsigned char *sig_msg, unsigned long long smlen, const unsigned char *pk) {

    if (smlen < CRYPTO_BYTES)
        goto badsig;

    *mlen = smlen - CRYPTO_BYTES;

    if (mirath_verify((uint8_t*)sig_msg + CRYPTO_BYTES, *mlen, (uint8_t*)sig_msg, (uint8_t*)pk) != EXIT_SUCCESS) {
        goto badsig;
    } else {
        /* All good, copy msg, return 0 */
        for (size_t i = 0; i < *mlen; ++i) msg[i] = sig_msg[CRYPTO_BYTES + i];
        return 0;
    }

    badsig:
    /* Signature verification failed */
    *mlen = -1;
    return -1;
}
