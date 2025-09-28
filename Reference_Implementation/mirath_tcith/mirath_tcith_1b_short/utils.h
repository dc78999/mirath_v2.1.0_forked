
#ifndef MIRATH_UTILS_H
#define MIRATH_UTILS_H

#include <stdint.h>
#include "framework_data_types.h"

/// taken from the kyber impl.
/// Description: Compare two arrays for equality in constant time.
///
/// Arguments:   const uint8_t *a: pointer to first byte array
///              const uint8_t *b: pointer to second byte array
///              size_t len:       length of the byte arrays
///
/// Returns 0 if the byte arrays are equal, 1 otherwise
static inline int hash_equal(const hash_t hash1, const hash_t hash2, const uint32_t len)
{
    uint8_t r = 0;
    uint32_t ret;

    for(uint32_t i=0; i < len; i++) {
        r |= hash1[i] ^ hash2[i];
    }

#ifdef __x86_64__
    ret = (-(uint64_t)r) >> 63;
#else
    ret = (-(uint32_t)r) >> 31;
#endif
    return ret == 0;
}

#endif //MIRATH_UTILS_H
