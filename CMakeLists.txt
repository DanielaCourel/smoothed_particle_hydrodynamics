cmake_minimum_required(VERSION 3.10)

project(sph)

# Use environment variables if set, otherwise find compilers
if(DEFINED ENV{CC})
    set(CMAKE_C_COMPILER $ENV{CC})
else()
    find_program(CMAKE_C_COMPILER NAMES gcc clang)
endif()

if(DEFINED ENV{CXX})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
else()
    find_program(CMAKE_CXX_COMPILER NAMES g++ clang++)
endif()

set(CMAKE_CXX_STANDARD 11)
# Modificar para agregar flags de compilación acá
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fPIC -Ofast -ffast-math -march=native -mavx2 -mfma")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O2 -ftree-vectorize -ffast-math -march=native -mfma -funroll -funroll-all-loops")

find_package(Qt5 COMPONENTS Core Gui OpenGL Widgets REQUIRED)

add_definitions(-D_USE_MATH_DEFINES)

set(OPENGL_LIBS GLU gomp GL)

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