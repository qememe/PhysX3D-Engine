cmake_minimum_required(VERSION 3.20)
project(PhysX3D VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags for optimization
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -Wall -Wextra -pedantic
        -march=native          # AMD Ryzen optimizations
        -mtune=native
        -O3                    # Maximum optimization
        -ffast-math
        -pthread
    )
endif()

# Library sources
add_library(physx3d STATIC
    src/math.cpp
    src/collision.cpp
    src/spatial.cpp
    src/core.cpp
    src/constraints.cpp
)

target_include_directories(physx3d PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Example executable
add_executable(physics_demo
    test/demo.cpp
)

target_link_libraries(physics_demo PRIVATE physx3d pthread)

# Install
install(TARGETS physx3d DESTINATION lib)
install(DIRECTORY include/physx3d DESTINATION include)
