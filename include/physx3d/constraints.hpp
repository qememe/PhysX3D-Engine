#pragma once

#include "math.hpp"

namespace physx3d {

class RigidBody;

// ============================================================================
// Base Constraint
// ============================================================================
class Constraint {
public:
    virtual ~Constraint() = default;
    virtual void solve(float dt) = 0;
    
protected:
    RigidBody* bodyA = nullptr;
    RigidBody* bodyB = nullptr;
};

// ============================================================================
// Distance Constraint (Spring)
// ============================================================================
class DistanceConstraint : public Constraint {
public:
    DistanceConstraint(RigidBody* a, RigidBody* b, float distance, float stiffness = 1.0f);
    void solve(float dt) override;
    
private:
    float targetDistance;
    float stiffness;
};

} // namespace physx3d
