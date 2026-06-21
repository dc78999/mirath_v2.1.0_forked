/**
 * @file rijndael.h
 * @brief AES-128 and Rijndael-256 implementation using ARM NEON Crypto Extensions
 */

#ifndef RIJNDAEL_H
#define RIJNDAEL_H

#include <arm_neon.h>
#include <string.h>

#define AES128_ROUNDS 10
#define RIJNDAEL256_ROUNDS 14

typedef struct {
    uint8x16_t keys[AES128_ROUNDS + 1];
} aes_128_round_keys_t;

typedef struct {
    uint8x16_t keys[(RIJNDAEL256_ROUNDS + 1)][2];
} rijndael_256_round_keys_t;

static inline uint8x16_t neon_slli_si128_4(uint8x16_t v) {
    return vextq_u8(vdupq_n_u8(0), v, 12);
}

static inline uint8x16_t _mm_shuffle_epi32_0xff(uint8x16_t v) {
    // 0xff = 11 11 11 11 -> broadcast word 3 to all words
    return vreinterpretq_u8_u32(vdupq_laneq_u32(vreinterpretq_u32_u8(v), 3));
}

static inline uint8x16_t _mm_shuffle_epi32_0xaa(uint8x16_t v) {
    // 0xaa = 10 10 10 10 -> broadcast word 2 to all words
    return vreinterpretq_u8_u32(vdupq_laneq_u32(vreinterpretq_u32_u8(v), 2));
}

static inline uint8x16_t aes_128_assist(uint8x16_t temp1, uint8x16_t temp2) {
    uint8x16_t temp3;
    temp3 = neon_slli_si128_4(temp1);
    temp1 = veorq_u8(temp1, temp3);
    temp3 = neon_slli_si128_4(temp3);
    temp1 = veorq_u8(temp1, temp3);
    temp3 = neon_slli_si128_4(temp3);
    temp1 = veorq_u8(temp1, temp3);
    temp1 = veorq_u8(temp1, temp2);
    return temp1;
}

static inline void aes_128_key_expansion(aes_128_round_keys_t *round_keys, const unsigned char *key) {
    uint8x16_t *Key_Schedule = round_keys->keys;
    
    // AVX2 _mm_set_epi32(0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d, 0x0c0f0e0d)
    const uint8_t shuffle_arr[16] = {
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c
    };
    uint8x16_t shuffle_mask = vld1q_u8(shuffle_arr);

    Key_Schedule[0] = vld1q_u8(key);

    uint8x16_t rcon = vreinterpretq_u8_u32(vdupq_n_u32(1));
    uint8x16_t tmp = vqtbl1q_u8(Key_Schedule[0], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[1] = aes_128_assist(Key_Schedule[0], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[2] = aes_128_assist(Key_Schedule[1], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[2], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[3] = aes_128_assist(Key_Schedule[2], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[3], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[4] = aes_128_assist(Key_Schedule[3], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[4], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[5] = aes_128_assist(Key_Schedule[4], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[5], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[6] = aes_128_assist(Key_Schedule[5], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[6], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[7] = aes_128_assist(Key_Schedule[6], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[7], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[8] = aes_128_assist(Key_Schedule[7], tmp);

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0x1b));
    tmp = vqtbl1q_u8(Key_Schedule[8], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[9] = aes_128_assist(Key_Schedule[8], tmp);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[9], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    Key_Schedule[10] = aes_128_assist(Key_Schedule[9], tmp);
}

static inline void aes_128_encrypt(unsigned char *out, const unsigned char *in,
                                   const aes_128_round_keys_t *Key_Schedule) {
    const uint8x16_t *KS = Key_Schedule->keys;
    uint8x16_t data = vld1q_u8(in);

    // Initial key addition. Since aeseq_u8 also adds the key, we don't need a separate XOR for round 0.
    // Wait, AVX2 does data = data ^ KS[0], then aesenc(data, KS[1]).
    // For ARM: aeseq_u8(data, KS[0]) does AddRoundKey(KS[0]) -> SubBytes -> ShiftRows
    // Then vaesmcq_u8 adds MixColumns.
    // So the sequence matches perfectly:
    // loop i=0..8:
    //    data = vaesmcq_u8(vaeseq_u8(data, KS[i]))
    // end
    // data = vaeseq_u8(data, KS[9])
    // data = veorq_u8(data, KS[10])

    for (int j = 0; j < AES128_ROUNDS - 1; j++) {
        data = vaesmcq_u8(vaeseq_u8(data, KS[j]));
    }
    data = vaeseq_u8(data, KS[AES128_ROUNDS - 1]);
    data = veorq_u8(data, KS[AES128_ROUNDS]);

    vst1q_u8(out, data);
}

static inline void rijndael_256_assist(const uint8x16_t round_key_in[2], uint8x16_t temp1, uint8x16_t round_key_out[2]) {
    uint8x16_t t1, t2, t3, t4;

    t1 = round_key_in[0];
    t3 = round_key_in[1];
    t2 = temp1;

    t2 = _mm_shuffle_epi32_0xff(t2);
    t4 = neon_slli_si128_4(t1);
    t1 = veorq_u8(t1, t4);
    t4 = neon_slli_si128_4(t4);
    t1 = veorq_u8(t1, t4);
    t4 = neon_slli_si128_4(t4);
    t1 = veorq_u8(t1, t4);
    t1 = veorq_u8(t1, t2);

    round_key_out[0] = t1;

    t4 = _mm_shuffle_epi32_0xff(t1);
    t4 = vaeseq_u8(t4, vdupq_n_u8(0));

    t2 = _mm_shuffle_epi32_0xaa(t4);
    t4 = neon_slli_si128_4(t3);
    t3 = veorq_u8(t3, t4);
    t4 = neon_slli_si128_4(t4);
    t3 = veorq_u8(t3, t4);
    t4 = neon_slli_si128_4(t4);
    t3 = veorq_u8(t3, t4);
    t3 = veorq_u8(t3, t2);

    round_key_out[1] = t3;
}

static inline void rijndael_256_key_expansion(rijndael_256_round_keys_t *round_keys, const unsigned char *key) {
    uint8x16_t(*Key_Schedule)[2] = round_keys->keys;
    
    const uint8_t shuffle_arr[16] = {
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c, 
        0x0d, 0x0e, 0x0f, 0x0c
    };
    uint8x16_t shuffle_mask = vld1q_u8(shuffle_arr);

    Key_Schedule[0][0] = vld1q_u8(key);
    Key_Schedule[0][1] = vld1q_u8(key + 16);

    uint8x16_t rcon = vreinterpretq_u8_u32(vdupq_n_u32(1));
    uint8x16_t tmp = vqtbl1q_u8(Key_Schedule[0][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[0], tmp, Key_Schedule[1]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[1][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[1], tmp, Key_Schedule[2]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[2][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[2], tmp, Key_Schedule[3]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[3][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[3], tmp, Key_Schedule[4]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[4][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[4], tmp, Key_Schedule[5]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[5][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[5], tmp, Key_Schedule[6]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[6][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[6], tmp, Key_Schedule[7]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[7][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[7], tmp, Key_Schedule[8]);

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0x1B));
    tmp = vqtbl1q_u8(Key_Schedule[8][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[8], tmp, Key_Schedule[9]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[9][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[9], tmp, Key_Schedule[10]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[10][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[10], tmp, Key_Schedule[11]);

    rcon = vreinterpretq_u8_u32(vshlq_n_u32(vreinterpretq_u32_u8(rcon), 1));
    tmp = vqtbl1q_u8(Key_Schedule[11][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[11], tmp, Key_Schedule[12]);

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0xAB));
    tmp = vqtbl1q_u8(Key_Schedule[12][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[12], tmp, Key_Schedule[13]);

    rcon = vreinterpretq_u8_u32(vdupq_n_u32(0x4D));
    tmp = vqtbl1q_u8(Key_Schedule[13][1], shuffle_mask);
    tmp = veorq_u8(vaeseq_u8(tmp, vdupq_n_u8(0)), rcon);
    rijndael_256_assist(Key_Schedule[13], tmp, Key_Schedule[14]);
}

static inline void rijndael_256_encrypt(unsigned char *out, const unsigned char *in, const rijndael_256_round_keys_t *Key_Schedule) {
    const uint8_t RIJNDAEL256_MASK_ARR[16] = {
        0x00, 0x01, 0x06, 0x07,
        0x04, 0x05, 0x0a, 0x0b,
        0x08, 0x09, 0x0e, 0x0f,
        0x0c, 0x0d, 0x02, 0x03
    };
    uint8x16_t RIJNDAEL256_MASK = vld1q_u8(RIJNDAEL256_MASK_ARR);
    
    const uint8_t BLEND_MASK_ARR[16] = {
        0x00, 0xFF, 0xFF, 0xFF,
        0x00, 0x00, 0xFF, 0xFF,
        0x00, 0x00, 0xFF, 0xFF,
        0x00, 0x00, 0x00, 0xFF
    };
    uint8x16_t BLEND_MASK = vld1q_u8(BLEND_MASK_ARR);

    const uint8x16_t(*KS)[2] = Key_Schedule->keys;
    uint8x16_t data1 = vld1q_u8(in);
    uint8x16_t data2 = vld1q_u8(in + 16);

    data1 = veorq_u8(data1, KS[0][0]);
    data2 = veorq_u8(data2, KS[0][1]);

    int j;
    for (j = 1; j < RIJNDAEL256_ROUNDS; j++) {
        uint8x16_t tmp1 = vbslq_u8(BLEND_MASK, data2, data1);
        uint8x16_t tmp2 = vbslq_u8(BLEND_MASK, data1, data2);
        
        tmp1 = vqtbl1q_u8(tmp1, RIJNDAEL256_MASK);
        tmp2 = vqtbl1q_u8(tmp2, RIJNDAEL256_MASK);
        
        data1 = vaesmcq_u8(vaeseq_u8(tmp1, vdupq_n_u8(0)));
        data1 = veorq_u8(data1, KS[j][0]);
        data2 = vaesmcq_u8(vaeseq_u8(tmp2, vdupq_n_u8(0)));
        data2 = veorq_u8(data2, KS[j][1]);
    }

    uint8x16_t tmp1 = vbslq_u8(BLEND_MASK, data2, data1);
    uint8x16_t tmp2 = vbslq_u8(BLEND_MASK, data1, data2);
    tmp1 = vqtbl1q_u8(tmp1, RIJNDAEL256_MASK);
    tmp2 = vqtbl1q_u8(tmp2, RIJNDAEL256_MASK);

    data1 = veorq_u8(vaeseq_u8(tmp1, vdupq_n_u8(0)), KS[j][0]);
    data2 = veorq_u8(vaeseq_u8(tmp2, vdupq_n_u8(0)), KS[j][1]);

    vst1q_u8(out, data1);
    vst1q_u8(out + 16, data2);
}

#endif //RIJNDAEL_H
