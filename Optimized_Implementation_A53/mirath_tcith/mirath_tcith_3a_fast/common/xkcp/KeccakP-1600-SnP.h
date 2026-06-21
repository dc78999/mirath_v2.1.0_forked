/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

Keccak-p permutation state definition for ARM (generic 64-bit).

When using KeccakP-1600-opt64.c, the state is simply 25 uint64_t lanes.
*/

#ifndef _KeccakP_1600_SnP_h_
#define _KeccakP_1600_SnP_h_

#include <stddef.h>
#include <stdint.h>

/* Generic 64-bit optimized state layout (25 lanes of 64 bits) */
typedef struct {
    ALIGN(8) uint64_t A[25];
} KeccakP1600_plain64_state;

typedef KeccakP1600_plain64_state KeccakP1600_state;

#define KeccakP1600_implementation      "ARMv8-A 64-bit optimized implementation"

#define KeccakP1600_StaticInitialize()
void KeccakP1600_Initialize(KeccakP1600_state *state);
void KeccakP1600_AddByte(KeccakP1600_state *state, unsigned char data, unsigned int offset);
void KeccakP1600_AddBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_OverwriteBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_OverwriteWithZeroes(KeccakP1600_state *state, unsigned int byteCount);
void KeccakP1600_Permute_Nrounds(KeccakP1600_state *state, unsigned int nrounds);
void KeccakP1600_Permute_12rounds(KeccakP1600_state *state);
void KeccakP1600_Permute_24rounds(KeccakP1600_state *state);
void KeccakP1600_ExtractBytes(const KeccakP1600_state *state, unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_ExtractAndAddBytes(const KeccakP1600_state *state, const unsigned char *input, unsigned char *output, unsigned int offset, unsigned int length);
size_t KeccakF1600_FastLoop_Absorb(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen);

#endif
