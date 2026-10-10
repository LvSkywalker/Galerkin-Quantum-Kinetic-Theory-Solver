#!/bin/sh
# Offline fallback when CMake is unavailable. CXX may select any C++17 compiler.
set -eu
cd "$(dirname "$0")"
mkdir -p build
SDK_FLAGS=
if [ -z "${CXX:-}" ] && [ -x /Library/Developer/CommandLineTools/usr/bin/clang++ ]; then
    CXX=/Library/Developer/CommandLineTools/usr/bin/clang++
    export DEVELOPER_DIR=/Library/Developer/CommandLineTools
    SDK_FLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
else
    CXX=${CXX:-c++}
fi
for file in src/numerics.cpp src/collisions.cpp src/solver.cpp src/dynamics.cpp src/thermodynamics.cpp src/fourier.cpp src/workflows.cpp src/figures.cpp src/main.cpp tests/tests.cpp; do
    object="build/$(basename "$file" .cpp).o"
    if [ ! -f "$object" ] || [ "$file" -nt "$object" ] || [ include/qtransport/core.hpp -nt "$object" ] || [ build.sh -nt "$object" ]; then
        echo "Compiling $file"
        "$CXX" $SDK_FLAGS -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Iinclude -isystem third_party/eigen -isystem third_party -c "$file" -o "$object"
    fi
done
"$CXX" $SDK_FLAGS build/numerics.o build/collisions.o build/solver.o build/dynamics.o build/thermodynamics.o build/fourier.o build/workflows.o build/figures.o build/main.o build/tests.o -o build/qtransport
echo "Built build/qtransport"
