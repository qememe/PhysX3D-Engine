#pragma once

#include "math.hpp"
#include "collision.hpp"
#include "spatial.hpp"
#include <vector>
#include <memory>
#include <thread>
#include <atomic>

namespace physx3d {

// ============================================================================
// RigidBody - Dynamic or static physics object
// ============================================================================
class RigidBody {
public:
    RigidBody() = default;
    
    // Getters
    const Vec3& getPosition() const { return position; }
    const Quaternion& getRotation() const { return rotation; }
    const Vec3& getVelocity() const { return velocity; }
    const Vec3& getAngularVelocity() const { return angularVelocity; }
    float getMass() const { return mass; }
    float getInvMass() const { return invMass; }
    bool isStatic() const { return isStaticBody; }
    bool isSleeping() const { return sleeping; }
    const Shape& getShape() const { return shape; }
    
    // Setters
    void setPosition(const Vec3& pos) { position = pos; }
    void setRotation(const Quaternion& rot) { rotation = rot.normalized(); }
    void setVelocity(const Vec3& vel) { velocity = vel; }
    void setAngularVelocity(const Vec3& angVel) { angularVelocity = angVel; }
    void setMass(float m);
    void setStatic(bool s) { isStaticBody = s; if (s) { invMass = 0; } }
    void setShape(const Shape& s) { shape = s; updateInertia(); }
    void setRestitution(float r) { restitution = r; }
    void setFriction(float f) { friction = f; }
    
    // Physics operations
    void applyForce(const Vec3& force) { if (!isStaticBody) this->force += force; }
    void applyImpulse(const Vec3& impulse) { if (!isStaticBody) velocity += impulse * invMass; }
    void applyTorque(const Vec3& torque) { if (!isStaticBody) this->torque += torque; }
    
    void applyImpulseAtPoint(const Vec3& impulse, const Vec3& point) {
        if (isStaticBody) return;
        velocity += impulse * invMass;
        Vec3 r = point - position;
        angularVelocity += invInertiaWorld * r.cross(impulse);
    }
    
    void integrate(float dt);
    void updateInertia();
    
    AABB getAABB() const {
        return computeAABB(shape, position, rotation);
    }
    
    friend class Constraint;
    friend class DistanceConstraint;

private:
    friend class PhysicsWorld;
    friend class DistanceConstraint;
    
    // Transform
    Vec3 position{0, 0, 0};
    Quaternion rotation;
    
    // Linear motion
    Vec3 velocity{0, 0, 0};
    Vec3 force{0, 0, 0};
    float mass = 1.0f;
    float invMass = 1.0f;
    
    // Angular motion
    Vec3 angularVelocity{0, 0, 0};
    Vec3 torque{0, 0, 0};
    Mat3 inertia = Mat3::identity();
    Mat3 invInertia = Mat3::identity();
    Mat3 invInertiaWorld = Mat3::identity();
    
    // Material properties
    float restitution = 0.5f; // Bounciness
    float friction = 0.5f;
    
    // Shape
    Shape shape = Sphere(1.0f);
    
    // State
    bool isStaticBody = false;
    bool sleeping = false;
    float sleepTimer = 0.0f;
};

// ============================================================================
// PhysicsWorld - Main simulation environment
// ============================================================================
class PhysicsWorld {
public:
    explicit PhysicsWorld(const Vec3& gravity = Vec3(0, -9.81f, 0));
    ~PhysicsWorld();
    
    // Body management
    RigidBody* createRigidBody();
    void removeRigidBody(RigidBody* body);
    
    // Simulation
    void step(float dt);
    
    // Configuration
    void setGravity(const Vec3& g) { gravity = g; }
    Vec3 getGravity() const { return gravity; }
    
    void setThreadCount(int count) { threadCount = std::max(1, count); }
    
    // Queries
    bool raycast(const Vec3& origin, const Vec3& direction, float maxDist, 
                 RigidBody*& hitBody, Vec3& hitPoint) const;
    
    friend class Constraint;
    friend class DistanceConstraint;

private:
    std::vector<std::unique_ptr<RigidBody>> bodies;
    Vec3 gravity;
    BVH bvh;
    
    int threadCount = std::thread::hardware_concurrency();
    
    // Simulation steps
    void applyGravity(float dt);
    void integrateVelocities(float dt);
    void broadPhase(std::vector<std::pair<int, int>>& pairs);
    void narrowPhase(const std::vector<std::pair<int, int>>& pairs, 
                     std::vector<Contact>& contacts);
    void solveContacts(std::vector<Contact>& contacts, float dt);
    void updateSleeping(float dt);
};

} // namespace physx3d
