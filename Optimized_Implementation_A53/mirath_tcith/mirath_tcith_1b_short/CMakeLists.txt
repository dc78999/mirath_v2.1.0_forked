cmake_minimum_required(VERSION 3.22)
project(mirath_tcith_1b_short C ASM)
set(CMAKE_C_STANDARD 11)

# define a default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

# A53 specific compiler flags
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -Wall -Wextra -march=armv8-a+crypto+crc -mtune=cortex-a53 -fomit-frame-pointer -funroll-loops -ftree-vectorize")

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

# We removed the OpenSSL requirement since rng.c uses /dev/urandom directly for cross-compilation

set(COMMON_SOURCES
    common/hash_sha3_128.c
    common/xkcp/KeccakHash.c
    common/xkcp/KeccakP-1600-opt64.c
)

set(MIRATH_SOURCES
    mirath_matrix_ff.c
    mirath_ggm_tree.c
    mirath_keygen.c
    tcith/mirath_tcith.c
    tcith/mirath_sign.c
    tcith/mirath_verify.c
    tcith/mirath_parsing.c
)

add_executable(PQCgenKAT_sign_${PROJECT_NAME} PQCgenKAT_sign.c sign.c rng.c ${MIRATH_SOURCES} ${COMMON_SOURCES})
target_include_directories(PQCgenKAT_sign_${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" common common/xkcp optimized/neon optimized/neon/arith tcith)

add_executable(mirath_tcith_1b_short main.c sign.c rng.c ${MIRATH_SOURCES} ${COMMON_SOURCES})
target_include_directories(mirath_tcith_1b_short PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" common common/xkcp optimized/neon optimized/neon/arith tcith)
