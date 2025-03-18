cmake_minimum_required(VERSION 3.10)

project(sph)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

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

target_include_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} Qt5::Core Qt5::Gui Qt5::OpenGL Qt5::Widgets ${OPENGL_LIBS})
