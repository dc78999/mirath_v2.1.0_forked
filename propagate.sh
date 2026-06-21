#!/bin/bash
set -e

src_dir="Optimized_Implementation/mirath_tcith"
dst_dir="Optimized_Implementation_A53/mirath_tcith"

variants=(
    "mirath_tcith_1a_fast" "mirath_tcith_1b_fast" "mirath_tcith_1a_short" "mirath_tcith_1b_short"
    "mirath_tcith_3a_fast" "mirath_tcith_3b_fast" "mirath_tcith_3a_short" "mirath_tcith_3b_short"
    "mirath_tcith_5a_fast" "mirath_tcith_5b_fast" "mirath_tcith_5a_short" "mirath_tcith_5b_short"
)

ref_dir="$dst_dir/mirath_tcith_1a_fast"

for var in "${variants[@]}"; do
    if [ "$var" == "mirath_tcith_1a_fast" ]; then continue; fi
    echo "Processing $var..."
    
    cp -r "$src_dir/$var" "$dst_dir/"
    
    cp -r "$ref_dir/optimized/neon" "$dst_dir/$var/optimized/"
    rm -rf "$dst_dir/$var/common"
    cp -r "$ref_dir/common" "$dst_dir/$var/"
    cp "$ref_dir/rng.c" "$dst_dir/$var/"
    cp "$ref_dir/main.c" "$dst_dir/$var/"
    cp "$ref_dir/mirath_arith.h" "$dst_dir/$var/"
    
    sed "s/mirath_tcith_1a_fast/$var/g" "$ref_dir/CMakeLists.txt" > "$dst_dir/$var/CMakeLists.txt"
done

# create a top-level CMakeLists.txt to build all 12 variants easily
cat << EOF > "Optimized_Implementation_A53/CMakeLists.txt"
cmake_minimum_required(VERSION 3.22)
project(mirath_tcith_all C ASM)

add_subdirectory(mirath_tcith/mirath_tcith_1a_fast)
add_subdirectory(mirath_tcith/mirath_tcith_1b_fast)
add_subdirectory(mirath_tcith/mirath_tcith_1a_short)
add_subdirectory(mirath_tcith/mirath_tcith_1b_short)
add_subdirectory(mirath_tcith/mirath_tcith_3a_fast)
add_subdirectory(mirath_tcith/mirath_tcith_3b_fast)
add_subdirectory(mirath_tcith/mirath_tcith_3a_short)
add_subdirectory(mirath_tcith/mirath_tcith_3b_short)
add_subdirectory(mirath_tcith/mirath_tcith_5a_fast)
add_subdirectory(mirath_tcith/mirath_tcith_5b_fast)
add_subdirectory(mirath_tcith/mirath_tcith_5a_short)
add_subdirectory(mirath_tcith/mirath_tcith_5b_short)
EOF

echo "Done propagating."
