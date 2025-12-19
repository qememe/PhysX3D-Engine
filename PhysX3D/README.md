# PhysX3D - High-Performance 3D Physics Engine

A modern C++20 physics engine optimized for multi-core CPUs (AMD Ryzen 12-thread tested).

## Features

### Core Physics
- **Rigid Body Dynamics**: Semi-implicit Euler integration
- **Collision Detection**: 
  - Broad-phase: BVH (Bounding Volume Hierarchy)
  - Narrow-phase: SAT, sphere-sphere, sphere-box, box-box
- **Collision Resolution**: Sequential impulse solver with friction and restitution
- **Shapes**: Sphere, Box, Capsule
- **Constraints**: Distance constraints (springs)

### Performance
- Multi-threaded collision detection (up to 12 threads)
- Cache-friendly data layouts
- SIMD-ready math library (AVX/SSE)
- Spatial partitioning (BVH) for O(log n) broad-phase
- Sleeping bodies optimization
- Target: >1000 rigid bodies @ 60 FPS

## Architecture

Structure:
- math.hpp: Vec3, Mat3, Quaternion (SIMD-aligned)
- collision.hpp: Shape types, AABB, collision detection
- spatial.hpp: BVH for spatial partitioning
- core.hpp: RigidBody, PhysicsWorld
- constraints.hpp: Distance constraint, joints

## Quick Start

See test/demo.cpp for example usage.

## Build

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j12
./physics_demo

### Requirements
- GCC 15.2+ or Clang 21+
- C++20 support
- CMake 3.20+

## Performance

System: AMD Ryzen 5 7500F (12 threads), 31GB RAM
Compiler: GCC 15.2.1 -O3 -march=native

Target: 1000+ rigid bodies at 60 FPS

## License

MIT License

## Author

Senior Backend Developer (C++/Python)
Worker ID: chat1
Built with: GCC 15.2.1, C++20, AMD Ryzen optimization
