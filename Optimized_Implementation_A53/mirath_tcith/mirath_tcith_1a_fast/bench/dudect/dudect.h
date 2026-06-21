#ifndef DUDECT_H
#define DUDECT_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>

typedef struct {
    double mean[2];
    double m2[2];
    uint64_t n[2];
} dudect_ctx_t;

static inline void dudect_init(dudect_ctx_t *ctx) {
    ctx->mean[0] = 0.0; ctx->mean[1] = 0.0;
    ctx->m2[0] = 0.0; ctx->m2[1] = 0.0;
    ctx->n[0] = 0; ctx->n[1] = 0;
}

static inline void dudect_push(dudect_ctx_t *ctx, uint64_t class_label, double value) {
    int c = class_label ? 1 : 0;
    ctx->n[c]++;
    double delta = value - ctx->mean[c];
    ctx->mean[c] += delta / ctx->n[c];
    double delta2 = value - ctx->mean[c];
    ctx->m2[c] += delta * delta2;
}

static inline double dudect_t_test(dudect_ctx_t *ctx) {
    if (ctx->n[0] < 2 || ctx->n[1] < 2) return 0.0;
    double var0 = ctx->m2[0] / (ctx->n[0] - 1);
    double var1 = ctx->m2[1] / (ctx->n[1] - 1);
    double num = ctx->mean[0] - ctx->mean[1];
    double den = sqrt(var0 / ctx->n[0] + var1 / ctx->n[1]);
    if (den == 0.0) return 0.0;
    return num / den;
}

static inline int dudect_report(dudect_ctx_t *ctx) {
    double t = dudect_t_test(ctx);
    double max_t = 4.5; // Roughly p=0.00001
    printf("dudect: t-statistic = %f ", t);
    if (fabs(t) > max_t) {
        printf("-> LEAK DETECTED (t > %f)\n", max_t);
        return 1;
    } else {
        printf("-> PASS\n");
        return 0;
    }
}

#endif // DUDECT_H
