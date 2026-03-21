cmake_minimum_required(VERSION 3.25)
project(planets VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for tooling (clangd, etc.)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler warnings
add_compile_options(
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# Library target (shared between app and tests)
add_library(planets_lib
    src/simulator.cpp
)
target_include_directories(planets_lib PUBLIC include)

# Main executable
add_executable(planets app/main.cpp)
target_link_libraries(planets PRIVATE planets_lib)

# Tests
enable_testing()
add_subdirectory(tests)
