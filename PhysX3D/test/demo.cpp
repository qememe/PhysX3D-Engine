#include <physx3d/core.hpp>
#include <iostream>
#include <chrono>

using namespace physx3d;

int main() {
    std::cout << "PhysX3D Demo - Falling Cubes\n";
    std::cout << "============================\n\n";
    
    // Create world
    PhysicsWorld world(Vec3(0, -9.81f, 0));
    
    // Create ground (static)
    auto ground = world.createRigidBody();
    ground->setShape(Box(Vec3(50, 1, 50)));
    ground->setPosition(Vec3(0, -1, 0));
    ground->setStatic(true);
    
    // Create falling cubes
    std::vector<RigidBody*> cubes;
    for (int i = 0; i < 5; ++i) {
        auto cube = world.createRigidBody();
        cube->setShape(Box(Vec3(0.5f, 0.5f, 0.5f)));
        cube->setPosition(Vec3(i * 1.5f - 3.0f, 10 + i * 2, 0));
        cube->setMass(1.0f);
        cube->setRestitution(0.3f);
        cubes.push_back(cube);
    }
    
    // Simulation loop
    const float dt = 1.0f / 60.0f;
    const int steps = 300; // 5 seconds
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    for (int step = 0; step < steps; ++step) {
        world.step(dt);
        
        // Print status every 30 steps (0.5s)
        if (step % 30 == 0) {
            std::cout << "Time: " << (step * dt) << "s\n";
            for (size_t i = 0; i < cubes.size(); ++i) {
                Vec3 pos = cubes[i]->getPosition();
                Vec3 vel = cubes[i]->getVelocity();
                std::cout << "  Cube " << i << ": pos=" << pos 
                         << " vel_y=" << vel.y << "\n";
            }
            std::cout << "\n";
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    std::cout << "\n============================\n";
    std::cout << "Simulation completed in " << duration.count() << "ms\n";
    std::cout << "Performance: " << (steps * 1000.0f / duration.count()) << " FPS\n";
    
    return 0;
}
