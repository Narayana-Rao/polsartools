cmake_minimum_required(VERSION 3.10)
project(polsartools)

set(CMAKE_CXX_STANDARD 14)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
add_library(polsartools_cpp STATIC
    src/processing.cpp
)

# Python bindings
add_library(polsartools_python_bindings MODULE
    python_bindings/bindings.cpp
)
target_link_libraries(polsartools_python_bindings polsartools_cpp)
