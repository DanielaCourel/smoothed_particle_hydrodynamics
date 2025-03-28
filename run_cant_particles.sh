#!/bin/bash

# Directorio de compilación
BUILD_DIR="build"
EXECUTABLE="sph"  # Nombre del ejecutable según CMakeLists.txt

# Limpiar y preparar el directorio de compilación
rm -rf $BUILD_DIR
mkdir $BUILD_DIR

CC="gcc"
CXX="g++"

echo "================================================"
echo "Compilando con: $CC / $CXX"
echo "================================================"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="outperf_${CC}_${TIMESTAMP}"

mkdir -p "$OUTPUT_DIR"

FLAGS="-O3 -march=native -mavx2"

echo "Compilando con FLAGS: $FLAGS"

rm -rf $BUILD_DIR/CMakeCache.txt $BUILD_DIR/CMakeFiles

# Configurar y compilar con el compilador seleccionado 
cmake -B $BUILD_DIR -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_CXX_FLAGS="$FLAGS"
cmake --build $BUILD_DIR

# Ejecutar y guardar salida
for i in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096;
do
    echo "Ejecutando con FLAGS: $FLAGS y $i particulas"
    perf stat -e cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses,L1-icache-loads,L1-icache-load-misses,dTLB-loads,dTLB-load-misses,iTLB-loads,iTLB-load-misses ./$BUILD_DIR/$EXECUTABLE r $i > "$OUTPUT_DIR/output_$(echo $FLAGS | tr ' ' '_')_$i.txt" 2>&1
done