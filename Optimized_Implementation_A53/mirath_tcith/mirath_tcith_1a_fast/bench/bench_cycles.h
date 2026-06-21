#ifndef BENCH_CYCLES_H
#define BENCH_CYCLES_H

#include <stdint.h>

static inline uint64_t cpucyclesStart(void) {
    uint64_t val;
    __asm__ __volatile__("mrs %0, cntvct_el0" : "=r" (val));
    return val;
}

static inline uint64_t cpucyclesStop(void) {
    uint64_t val;
    __asm__ __volatile__("mrs %0, cntvct_el0" : "=r" (val));
    return val;
}

#endif // BENCH_CYCLES_H
