cmake_minimum_required(VERSION 3.0)

project(dynamic_color)

find_package(OpenMP REQUIRED)

set (USE_PYBIND11 TRUE)
# set (USE_PYBIND11 FALSE)

if(USE_PYBIND11)
    set(SOURCES
        src/dynamic_color.cpp
        src/evaluator_discri.cpp
        src/palette_evaluator.cpp
        src/grid_data.cpp
        src/scatter_data.cpp
        src/line_data.cpp
        src/generator.cpp)
else()
    set(SOURCES
        src/main.cpp
        src/evaluator_discri.cpp
        src/palette_evaluator.cpp
        src/grid_data.cpp
        src/scatter_data.cpp
        src/line_data.cpp
        src/generator.cpp)
endif()


# Add static files
configure_file(static/c3_data.json ${CMAKE_CURRENT_BINARY_DIR}/static/c3_data.json COPYONLY)

# Add header files
include_directories(include)

# set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")

set(CMAKE_BUILD_TYPE "Release")
# use optimization level 3
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


if(USE_PYBIND11)
    add_subdirectory(pybind11)
    pybind11_add_module(dynamic_color MODULE ${SOURCES})
    target_compile_options(dynamic_color PUBLIC -fopenmp)
    # remove -pg
    # set_target_properties(dynamic_color PROPERTIES LINK_FLAGS -pg)
    target_compile_features(dynamic_color PUBLIC cxx_std_17)
    target_link_libraries(dynamic_color PRIVATE OpenMP::OpenMP_CXX)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -fopenmp")
    add_executable(dynamic_color ${SOURCES})
    target_compile_features(dynamic_color PUBLIC cxx_std_17)
    target_link_libraries(dynamic_color PRIVATE OpenMP::OpenMP_CXX)
endif()

