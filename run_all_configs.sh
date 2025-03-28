#!/bin/bash

# Directorio de compilación
BUILD_DIR="build"
EXECUTABLE="sph"  # Nombre del ejecutable según CMakeLists.txt

# Lista de compiladores a probar
COMPILERS=(
    "gcc g++"
    "clang clang++"
)

# Diferentes configuraciones de compilación
CONFIGURATIONS=(
    "-O0 -march=native"
    "-O1 -march=native"
    "-O2 -march=native -mavx2"
    "-O3 -march=native -mavx2"
)

# Limpiar y preparar el directorio de compilación
rm -rf $BUILD_DIR
mkdir $BUILD_DIR

# Iterar sobre los compiladores
for COMPILER in "${COMPILERS[@]}"; do
    CC=$(echo $COMPILER | cut -d ' ' -f1)
    CXX=$(echo $COMPILER | cut -d ' ' -f2)

    echo "================================================"
    echo "Compilando con: $CC / $CXX"
    echo "================================================"

    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    OUTPUT_DIR="out_${CC}_${TIMESTAMP}"
    mkdir -p "$OUTPUT_DIR"

    for FLAGS in "${CONFIGURATIONS[@]}"; do
        echo "Compilando con FLAGS: $FLAGS"

        rm -rf $BUILD_DIR/CMakeCache.txt $BUILD_DIR/CMakeFiles

        # Configurar y compilar con el compilador seleccionado
        cmake -B $BUILD_DIR -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_CXX_FLAGS="$FLAGS"
        cmake --build $BUILD_DIR

        # Ejecutar y guardar salida
        echo "Ejecutando con FLAGS: $FLAGS"
        perf stat ./$BUILD_DIR/$EXECUTABLE r > "$OUTPUT_DIR/output_$(echo $FLAGS | tr ' ' '_').txt" 2>&1
    done

    echo "Los resultados fueron guardados en $OUTPUT_DIR"
done
