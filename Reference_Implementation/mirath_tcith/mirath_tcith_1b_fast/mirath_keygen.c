#include <stdint.h>

#include "mirath_matrix_ff.h"
#include "random.h"
#include "mirath_parsing.h"

int mirath_keygen(uint8_t *pk, uint8_t *sk) {
    seed_t seed_sk;
    seed_t seed_pk;
    ff_t S[MIRATH_VAR_FF_S_BYTES];
    ff_t C[MIRATH_VAR_FF_C_BYTES];
    ff_t H[MIRATH_VAR_FF_H_BYTES];
    ff_t y[MIRATH_VAR_FF_Y_BYTES];

    // step 1
    randombytes(seed_sk, MIRATH_SECURITY_BYTES);

    //step 2
    randombytes(seed_pk, MIRATH_SECURITY_BYTES);

    // step 3
    mirath_matrix_expand_seed_secret_matrix(S, C, seed_sk);

    // step 4
    mirath_matrix_expand_seed_public_matrix(H, seed_pk);

    // step 5
    mirath_matrix_compute_y(y, S, C, H);

    // step 6
    unparse_public_key(pk, seed_pk, y);

    // step 7
    unparse_secret_key(sk, seed_sk, seed_pk);

#if defined(_VERBOSE_)
    // ADD somethning here
    printf("\n+++++++ KEYGEN +++++++\n");
#endif

    return 0;
}
