# Mirath - Cortex-A53 Optimized Implementation

This directory contains the optimized implementation of **Mirath** (v2.1.0) targeting the **ARM Cortex-A53** architecture (e.g., Raspberry Pi Zero 2 W). 

## Optimizations

The reference and generic AVX2 implementations have been migrated and heavily optimized for the ARMv8-A architecture:

1. **Galois Field Arithmetic (NEON)**
   - All `__m256i` AVX2 operations have been ported to 128-bit ARM NEON `uint8x16_t` vectors.
   - Multiplication in GF(2⁴) and GF(2⁸) utilizes hardware-accelerated table lookups (`vqtbl1q_u8`) and fast bitwise XOR cascades (`veorq_u8`).
   - Constant-time bitwise selection leverages `vbslq_u8`.

2. **AES/Rijndael Key Expansion and Encryption (ARM Crypto Extensions)**
   - The tree expansion and PRG (Pseudo-Random Generator) logic utilizes the ARMv8 Cryptography Extensions.
   - Replaced software/lookup-based AES with hardware-accelerated, constant-time `vaeseq_u8` and `vaesmcq_u8` instructions.

3. **Keccak/SHA-3**
   - Integrates the 64-bit optimized Keccak permutation (`KeccakP-1600-opt64`) from XKCP, ensuring compatibility with the sequential x4 hashing APIs used in the TCitH tree generation.

4. **Constant-Time Guarantees**
   - The implementation strips out data-dependent branches and memory lookups. The GF arithmetic relies solely on constant-time NEON instructions.
   - Tested for constant-time vulnerabilities using a bundled `dudect` statistical evaluation harness (Welch's t-test).

## Prerequisites

To cross-compile this implementation from an x86_64 host to the AArch64 target, you will need the `aarch64-linux-gnu` toolchain installed:

```bash
# On Debian/Ubuntu
sudo apt-get install gcc-aarch64-linux-gnu g++-aarch64-linux-gnu
```

## Building

The build system uses CMake and a provided cross-compilation toolchain file. A top-level build will compile all 12 Mirath parameter sets simultaneously.

```bash
mkdir build && cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../cross-aarch64.cmake ..
make -j$(nproc)
```

This will generate the binaries (e.g., `mirath_tcith_1a_fast`, `PQCgenKAT_sign_mirath_tcith_1a_fast`) for all parameter variants.

## Benchmarking & Evaluation

Each parameter variant includes a dedicated `bench` module for performance and side-channel evaluation. These binaries must be run directly on the target ARM hardware (e.g., Raspberry Pi Zero 2 W).

### Cycle Benchmarking

The timing framework uses the native ARMv8 Generic Timer (`cntvct_el0` via `mrs` inline assembly) to accurately count cycles without requiring a patched kernel or `perf` access.

Run the main executable for your chosen variant to see median cycle counts for KeyGen, Sign, and Verify over 100 iterations:

```bash
./mirath_tcith_1a_fast
```

### Dudect (Side-Channel Evaluation)

To mathematically verify the constant-time execution of the `crypto_sign` operation, we've bundled a `dudect` harness that uses Welford's online algorithm to compute the Welch's t-test between fixed-message and random-message classes.

*(Note: You will need to compile the `ct_test.c` harness manually or add it to the CMake targets if you wish to run the empirical test).*

## Repository Structure

- `cross-aarch64.cmake`: The CMake cross-compilation toolchain configuration.
- `mirath_tcith/`: Contains the 12 parameter sets (`1a_fast`, `1a_short`, `3b_fast`, etc.).
  - `optimized/neon/`: Contains the core ARM NEON and Crypto extensions ports.
  - `bench/`: Contains the Dudect side-channel and benchmarking harnesses.
  - `tcith/`: The TCitH threshold signature core logic.
  - `common/`: XKCP Keccak and SHA3 wrappers.

## License
Refer to the root Mirath repository for licensing information.
