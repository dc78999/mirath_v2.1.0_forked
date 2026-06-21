#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>

static inline uint32_t aes_sub_word(uint32_t w) {
    uint32x4_t v = vdupq_n_u32(w);
    uint8x16_t v8 = vreinterpretq_u8_u32(v);
    v8 = vaeseq_u8(v8, vdupq_n_u8(0));
    return vgetq_lane_u32(vreinterpretq_u32_u8(v8), 0);
}

int main() {
    uint32_t w = 0x01020304;
    printf("w: %08x, sub_word: %08x\n", w, aes_sub_word(w));
    // expected: 
    // 01 -> 0x7c
    // 02 -> 0x77
    // 03 -> 0x7b
    // 04 -> 0xf2
    // If little endian, 0x01 is the lowest byte (b0). So 0x01020304 means b0=04, b1=03, b2=02, b3=01.
    // wait, w is uint32_t. In memory, it's 04 03 02 01 on little endian.
    // So 04 -> f2, 03 -> 7b, 02 -> 77, 01 -> 7c.
    // So output should be 0x7c777bf2.
    return 0;
}
