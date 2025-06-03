cmake_minimum_required(VERSION 3.10)

project(sph LANGUAGES CXX CUDA)

# -------------------------------------
# Selección del compilador
# -------------------------------------
if(DEFINED ENV{CC})
    set(CMAKE_C_COMPILER $ENV{CC})
endif()

if(DEFINED ENV{CXX})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
endif()

# Si aún no está definido, buscar los disponibles
if(NOT CMAKE_C_COMPILER)
    find_program(CMAKE_C_COMPILER NAMES icc icx gcc clang)
endif()

if(NOT CMAKE_CXX_COMPILER)
    find_program(CMAKE_CXX_COMPILER NAMES icpc icpx g++ clang++)
endif()

# -------------------------------------
# Estándar C++ y configuración general
# -------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -------------------------------------
# Flags por compilador
# -------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(STATUS "Usando compilador Intel: ${CMAKE_CXX_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -O3 -march=native")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVIDIA")
    message(STATUS "Usando compilador NVIDIA HPC: ${CMAKE_CXX_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -mp")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Usando compilador GCC: ${CMAKE_CXX_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -ffast-math -march=native -funroll-loops -fPIC")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Usando compilador Clang: ${CMAKE_CXX_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -ffast-math -march=native -fPIC")

else()
    message(WARNING "Compilador desconocido: ${CMAKE_CXX_COMPILER}. Se usarán flags por defecto.")
endif()

# -------------------------------------
# CUDA Configuration
# -------------------------------------
# Find CUDA toolkit
# find_package(CUDA REQUIRED) -> already done it on "LANGUAGE"...

# Set CUDA standard (e.g., C++11, C++14, C++17). This helps nvcc compile host code with the specified C++ standard.
set(CUDA_PROPAGATE_HOST_FLAGS OFF) # This prevents host compiler flags from being automatically propagated to nvcc.
set(CMAKE_CUDA_STANDARD 11) # Or whatever C++ standard you're using for your host code

# Add common CUDA compilation flags
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -gencode arch=sm_60")  # OLD
set(CMAKE_CUDA_ARCHITECTURES "61") # For "Pascal"
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")

# -------------------------------------
# Qt y OpenGL
# -------------------------------------
find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets REQUIRED)
add_definitions(-D_USE_MATH_DEFINES)

set(OPENGL_LIBS GLU gomp GL)

# -------------------------------------
# Archivos fuente
# -------------------------------------
set(SOURCES
    src/main.cpp
    src/sph.cpp
    src/particle.cpp
    src/visualization.cpp
    src/widget.cpp
    src/sphconfig.cpp
    src/vec3.cpp
)

set(HEADERS
    src/sph.h
    src/particle.h
    src/visualization.h
    src/widget.h
    src/sphconfig.h
    src/vec3.h
    src/cuda_utils.h        # Add CUDA header
)

set(FORMS
    src/widget.ui
)

# Add your CUDA source files here (e.g., .cu files)
set(CUDA_SOURCES
    src/cuda_utils.cu       # Dedicated CUDA source file
)

qt5_wrap_ui(UI_HEADERS ${FORMS})
qt5_wrap_cpp(MOC_SOURCES ${HEADERS})

add_executable(${PROJECT_NAME} ${SOURCES} ${MOC_SOURCES} ${UI_HEADERS} ${CUDA_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE src ${CMAKE_BINARY_DIR} ${CUDA_INCLUDE_DIRS}) # Add CUDA include directories
target_link_libraries(${PROJECT_NAME} Qt5::Core Qt5::Gui Qt5::OpenGL Qt5::Widgets ${OPENGL_LIBS} ${CUDA_LIBRARIES}) # Link CUDA libraries
