cmake_minimum_required(VERSION 3.10)

project(sph)

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
)

set(FORMS
    src/widget.ui
)

qt5_wrap_ui(UI_HEADERS ${FORMS})
qt5_wrap_cpp(MOC_SOURCES ${HEADERS})

add_executable(${PROJECT_NAME} ${SOURCES} ${MOC_SOURCES} ${UI_HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE src ${CMAKE_BINARY_DIR})
target_link_libraries(${PROJECT_NAME} Qt5::Core Qt5::Gui Qt5::OpenGL Qt5::Widgets ${OPENGL_LIBS})