#include "../include/physx3d/core.hpp"
#include <algorithm>
#include <cmath>

namespace physx3d {

// ============================================================================
// RigidBody Implementation
// ============================================================================

void RigidBody::setMass(float m) {
    mass = m;
    invMass = (m > 0 && !isStaticBody) ? 1.0f / m : 0.0f;
    updateInertia();
}

void RigidBody::updateInertia() {
    if (isStaticBody || mass <= 0) {
        invInertia = Mat3::scale(0);
        invInertiaWorld = Mat3::scale(0);
        return;
    }
    
    // Compute inertia tensor based on shape
    std::visit([&](auto&& s) {
        using T = std::decay_t<decltype(s)>;
        
        if constexpr (std::is_same_v<T, Sphere>) {
            float i = 0.4f * mass * s.radius * s.radius;
            inertia = Mat3::scale(i);
        }
        else if constexpr (std::is_same_v<T, Box>) {
            Vec3 he = s.halfExtents;
            float ix = (1.0f / 12.0f) * mass * (he.y * he.y + he.z * he.z) * 4;
            float iy = (1.0f / 12.0f) * mass * (he.x * he.x + he.z * he.z) * 4;
            float iz = (1.0f / 12.0f) * mass * (he.x * he.x + he.y * he.y) * 4;
            inertia = Mat3::diagonal(ix, iy, iz);
        }
        else if constexpr (std::is_same_v<T, Capsule>) {
            float i = 0.4f * mass * s.radius * s.radius; // Approximate
            inertia = Mat3::scale(i);
        }
    }, shape);
    
    invInertia = inertia.inverse();
    
    // Update world-space inverse inertia
    Mat3 rot = rotation.toMat3();
    invInertiaWorld = rot * invInertia * rot.transpose();
}

void RigidBody::integrate(float dt) {
    if (isStaticBody || sleeping) return;
    
    // Linear integration (semi-implicit Euler)
    velocity += force * invMass * dt;
    position += velocity * dt;
    
    // Angular integration
    angularVelocity += invInertiaWorld * torque * dt;
    
    // Update rotation
    Quaternion spin(0, angularVelocity.x, angularVelocity.y, angularVelocity.z);
    rotation += (spin * rotation) * (0.5f * dt);
    rotation.normalize();
    
    // Update world-space inertia
    updateInertia();
    
    // Clear forces
    force = Vec3(0, 0, 0);
    torque = Vec3(0, 0, 0);
    
    // Sleeping check
    float speedSq = velocity.lengthSq() + angularVelocity.lengthSq();
    if (speedSq < 0.01f) {
        sleepTimer += dt;
        if (sleepTimer > 1.0f) {
            sleeping = true;
        }
    } else {
        sleepTimer = 0;
    }
}

// ============================================================================
// PhysicsWorld Implementation
// ============================================================================

PhysicsWorld::PhysicsWorld(const Vec3& gravity) 
    : gravity(gravity) 
{
}

PhysicsWorld::~PhysicsWorld() = default;

RigidBody* PhysicsWorld::createRigidBody() {
    bodies.emplace_back(std::make_unique<RigidBody>());
    return bodies.back().get();
}

void PhysicsWorld::removeRigidBody(RigidBody* body) {
    bodies.erase(
        std::remove_if(bodies.begin(), bodies.end(),
            [body](const auto& ptr) { return ptr.get() == body; }),
        bodies.end()
    );
}

void PhysicsWorld::step(float dt) {
    // 1. Apply gravity
    applyGravity(dt);
    
    // 2. Broad-phase collision detection
    std::vector<std::pair<int, int>> pairs;
    broadPhase(pairs);
    
    // 3. Narrow-phase collision detection
    std::vector<Contact> contacts;
    narrowPhase(pairs, contacts);
    
    // 4. Solve contacts (iterative)
    for (int iter = 0; iter < 10; ++iter) {
        solveContacts(contacts, dt);
    }
    
    // 5. Integrate velocities
    integrateVelocities(dt);
    
    // 6. Update sleeping states
    updateSleeping(dt);
}

void PhysicsWorld::applyGravity(float dt) {
    for (auto& body : bodies) {
        if (!body->isStatic() && !body->isSleeping()) {
            body->applyForce(gravity * body->getMass());
        }
    }
}

void PhysicsWorld::integrateVelocities(float dt) {
    for (auto& body : bodies) {
        body->integrate(dt);
    }
}

void PhysicsWorld::broadPhase(std::vector<std::pair<int, int>>& pairs) {
    // Build BVH
    std::vector<AABB> aabbs;
    std::vector<int> indices;
    
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (!bodies[i]->isSleeping()) {
            aabbs.push_back(bodies[i]->getAABB());
            indices.push_back(i);
        }
    }
    
    if (!aabbs.empty()) {
        bvh.build(aabbs, indices);
        bvh.queryOverlaps(pairs);
    }
}

void PhysicsWorld::narrowPhase(const std::vector<std::pair<int, int>>& pairs, 
                               std::vector<Contact>& contacts)
{
    contacts.clear();
    
    for (const auto& [idxA, idxB] : pairs) {
        if (idxA >= bodies.size() || idxB >= bodies.size()) continue;
        
        RigidBody* a = bodies[idxA].get();
        RigidBody* b = bodies[idxB].get();
        
        // Skip if both static
        if (a->isStatic() && b->isStatic()) continue;
        
        Contact contact;
        if (detectCollision(
            a->getShape(), a->getPosition(), a->getRotation(),
            b->getShape(), b->getPosition(), b->getRotation(),
            contact))
        {
            contact.bodyA = idxA;
            contact.bodyB = idxB;
            contacts.push_back(contact);
            
            // Wake up sleeping bodies
            if (a->isSleeping()) a->sleeping = false;
            if (b->isSleeping()) b->sleeping = false;
        }
    }
}

void PhysicsWorld::solveContacts(std::vector<Contact>& contacts, float dt) {
    for (auto& contact : contacts) {
        if (contact.bodyA < 0 || contact.bodyB < 0) continue;
        if (contact.bodyA >= bodies.size() || contact.bodyB >= bodies.size()) continue;
        
        RigidBody* a = bodies[contact.bodyA].get();
        RigidBody* b = bodies[contact.bodyB].get();
        
        Vec3 rv = b->getVelocity() - a->getVelocity();
        float normalVel = rv.dot(contact.normal);
        
        // Don't resolve if separating
        if (normalVel > 0) continue;
        
        // Compute impulse
        float restitution = std::min(a->restitution, b->restitution);
        float j = -(1 + restitution) * normalVel;
        j /= (a->invMass + b->invMass);
        
        Vec3 impulse = contact.normal * j;
        
        // Apply impulse
        a->applyImpulse(-impulse);
        b->applyImpulse(impulse);
        
        // Position correction (Baumgarte stabilization)
        const float percent = 0.2f;
        const float slop = 0.01f;
        float correction = std::max(contact.penetration - slop, 0.0f) * percent;
        Vec3 correctionVec = contact.normal * correction;
        
        if (!a->isStatic()) {
            a->position -= correctionVec * (a->invMass / (a->invMass + b->invMass));
        }
        if (!b->isStatic()) {
            b->position += correctionVec * (b->invMass / (a->invMass + b->invMass));
        }
    }
}

void PhysicsWorld::updateSleeping(float dt) {
    // Already handled in RigidBody::integrate
}

bool PhysicsWorld::raycast(const Vec3& origin, const Vec3& direction, float maxDist, 
                           RigidBody*& hitBody, Vec3& hitPoint) const
{
    // TODO: Implement raycast
    return false;
}

} // namespace physx3d
