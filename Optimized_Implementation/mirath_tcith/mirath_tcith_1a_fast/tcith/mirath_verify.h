#ifndef MIRATH_TCITH_MIRATH_VERIFY_H
#define MIRATH_TCITH_MIRATH_VERIFY_H

#include <stddef.h>
#include <stdint.h>

int mirath_verify(uint8_t *msg, size_t msg_len, uint8_t *sig_msg, uint8_t *pk);

#endif //MIRATH_TCITH_MIRATH_VERIFY_H
