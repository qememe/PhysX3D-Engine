#include "../include/physx3d/constraints.hpp"
#include "../include/physx3d/core.hpp"

namespace physx3d {

DistanceConstraint::DistanceConstraint(RigidBody* a, RigidBody* b, float distance, float stiffness)
    : targetDistance(distance), stiffness(stiffness)
{
    bodyA = a;
    bodyB = b;
}

void DistanceConstraint::solve(float dt) {
    if (!bodyA || !bodyB) return;
    
    Vec3 delta = bodyB->getPosition() - bodyA->getPosition();
    float currentDist = delta.length();
    
    if (currentDist < 1e-6f) return;
    
    float error = currentDist - targetDistance;
    Vec3 correction = (delta / currentDist) * error * stiffness;
    
    float totalInvMass = bodyA->invMass + bodyB->invMass;
    if (totalInvMass < 1e-6f) return;
    
    if (!bodyA->isStatic()) {
        bodyA->position += correction * (bodyA->invMass / totalInvMass);
    }
    if (!bodyB->isStatic()) {
        bodyB->position -= correction * (bodyB->invMass / totalInvMass);
    }
}

} // namespace physx3d
