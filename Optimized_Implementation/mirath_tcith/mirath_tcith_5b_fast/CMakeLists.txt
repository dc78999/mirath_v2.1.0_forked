cmake_minimum_required(VERSION 3.22)
project(mirath_tcith_5b_fast C ASM)
set(CMAKE_C_STANDARD 11)

# define a default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -ftree-vectorize -funroll-loops -fomit-frame-pointer -fno-stack-protector -Wextra -mavx -mavx2 -mbmi -mbmi2 -maes")

if(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL x86_64)
    set(CMAKE_ASM_FLAGS "-x assembler-with-cpp -Wa,-defsym,old_gas_syntax=1 -Wa,-defsym,no_plt=1")
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if (EXISTS /usr/local/opt/openssl@3)
        set(OPENSSL_ROOT_DIR /usr/local/opt/openssl@3)
    elseif (EXISTS /usr/local/opt/openssl/)
        set(OPENSSL_ROOT_DIR /usr/local/opt/openssl/)
    endif()
endif()

# commit based on sha3 or rijndael
# by default, rijndael is used
string(TOLOWER "${OPT_COMMIT}" COMMIT_LOWER)
if (DEFINED OPT_COMMIT AND COMMIT_LOWER MATCHES sha3)
    add_compile_definitions(_SHA3_)
    message(STATUS "Commit based on: SHA3")
else ()
    message(STATUS "Commit based on: Rijndael")
endif ()

find_package(OpenSSL QUIET)

if (OPENSSL_FOUND)
    message(STATUS "OpenSSL: Version ${OPENSSL_VERSION}")
    message(STATUS "OpenSSL: include dir at ${OPENSSL_INCLUDE_DIR}")
    message(STATUS "OpenSSL: libraries at ${OPENSSL_LIBRARIES} ${OPENSSL_SSL_LIBRARIES}")
else ()
    message(FATAL_ERROR "OpenSSL NO FOUND")
endif()

set(COMMON_SOURCES
    common/hash_sha3_256.c;
    common/xkcp/KeccakHash.c;
    common/xkcp/KeccakP-1600-AVX2.s
    common/xkcp/avx2/KeccakHashtimes4.c
    common/xkcp/avx2/KeccakP-1600-times4-SIMD256.c
    common/xkcp/avx2/KeccakSpongetimes4.c

)

set(MIRATH_SOURCES
    mirath_matrix_ff.c
    mirath_ggm_tree.c
    mirath_keygen.c
    tcith/mirath_tcith.c
    tcith/mirath_sign.c
    tcith/mirath_verify.c
    tcith/mirath_parsing.c;
)

add_executable(PQCgenKAT_sign PQCgenKAT_sign.c sign.c rng.c ${MIRATH_SOURCES} ${COMMON_SOURCES})
target_include_directories(PQCgenKAT_sign PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" common common/xkcp common/xkcp/avx2 optimized/avx2 optimized/avx2/arith tcith)
target_link_libraries(PQCgenKAT_sign PRIVATE OpenSSL::Crypto)

add_executable(mirath_tcith_5b_fast main.c sign.c rng.c ${MIRATH_SOURCES} ${COMMON_SOURCES})
target_include_directories(mirath_tcith_5b_fast PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" common common/xkcp common/xkcp/avx2 optimized/avx2 optimized/avx2/arith tcith)
target_link_libraries(mirath_tcith_5b_fast PRIVATE OpenSSL::Crypto)
