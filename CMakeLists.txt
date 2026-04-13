cmake_minimum_required(VERSION 3.25)
project(planets VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_COLOR_DIAGNOSTICS ON)

include(FetchContent)
FetchContent_Declare(
    mp-units
    GIT_REPOSITORY https://github.com/mpusz/mp-units.git
    GIT_TAG        93d5ff2226e3ba26047bf12cb44c4620b8121dfe
    SOURCE_SUBDIR  src
)
set(MP_UNITS_API_CONTRACTS "NONE" CACHE STRING "" FORCE)
FetchContent_MakeAvailable(mp-units)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG        5.0
    GIT_SHALLOW    TRUE
)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_GAMES    OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(raylib)

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
target_link_libraries(planets_lib PUBLIC mp-units::mp-units)

# Main executable
add_executable(planets app/main.cpp app/visualization.cpp app/trails.cpp)
target_include_directories(planets PRIVATE app)
target_link_libraries(planets PRIVATE planets_lib raylib)

# Tests
enable_testing()
add_subdirectory(tests)
