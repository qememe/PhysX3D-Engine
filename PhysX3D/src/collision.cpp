#include "../include/physx3d/collision.hpp"
#include <algorithm>
#include <cmath>

namespace physx3d {

AABB computeAABB(const Shape& shape, const Vec3& position, const Quaternion& rotation) {
    return std::visit([&](auto&& s) -> AABB {
        using T = std::decay_t<decltype(s)>;
        
        if constexpr (std::is_same_v<T, Sphere>) {
            Vec3 r(s.radius, s.radius, s.radius);
            return AABB(position - r, position + r);
        }
        else if constexpr (std::is_same_v<T, Box>) {
            // Rotate box corners and find AABB
            Mat3 rot = rotation.toMat3();
            Vec3 corners[8];
            Vec3 he = s.halfExtents;
            corners[0] = rot * Vec3(-he.x, -he.y, -he.z);
            corners[1] = rot * Vec3( he.x, -he.y, -he.z);
            corners[2] = rot * Vec3(-he.x,  he.y, -he.z);
            corners[3] = rot * Vec3( he.x,  he.y, -he.z);
            corners[4] = rot * Vec3(-he.x, -he.y,  he.z);
            corners[5] = rot * Vec3( he.x, -he.y,  he.z);
            corners[6] = rot * Vec3(-he.x,  he.y,  he.z);
            corners[7] = rot * Vec3( he.x,  he.y,  he.z);
            
            AABB aabb;
            aabb.min = aabb.max = position + corners[0];
            for (int i = 1; i < 8; ++i) {
                aabb.expand(position + corners[i]);
            }
            return aabb;
        }
        else if constexpr (std::is_same_v<T, Capsule>) {
            float r = s.radius;
            float h = s.halfHeight;
            Vec3 extent(r, r + h, r);
            return AABB(position - extent, position + extent);
        }
        
        return AABB();
    }, shape);
}

bool detectCollision(
    const Shape& shapeA, const Vec3& posA, const Quaternion& rotA,
    const Shape& shapeB, const Vec3& posB, const Quaternion& rotB,
    Contact& contact)
{
    // Dispatch to specific collision functions
    if (std::holds_alternative<Sphere>(shapeA) && std::holds_alternative<Sphere>(shapeB)) {
        return sphereVsSphere(std::get<Sphere>(shapeA), posA,
                             std::get<Sphere>(shapeB), posB, contact);
    }
    else if (std::holds_alternative<Sphere>(shapeA) && std::holds_alternative<Box>(shapeB)) {
        return sphereVsBox(std::get<Sphere>(shapeA), posA, rotA,
                          std::get<Box>(shapeB), posB, rotB, contact);
    }
    else if (std::holds_alternative<Box>(shapeA) && std::holds_alternative<Sphere>(shapeB)) {
        bool hit = sphereVsBox(std::get<Sphere>(shapeB), posB, rotB,
                              std::get<Box>(shapeA), posA, rotA, contact);
        if (hit) {
            contact.normal = -contact.normal;
            std::swap(contact.bodyA, contact.bodyB);
        }
        return hit;
    }
    else if (std::holds_alternative<Box>(shapeA) && std::holds_alternative<Box>(shapeB)) {
        return boxVsBox(std::get<Box>(shapeA), posA, rotA,
                       std::get<Box>(shapeB), posB, rotB, contact);
    }
    
    return false; // Unsupported combination
}

bool sphereVsSphere(const Sphere& a, const Vec3& posA,
                    const Sphere& b, const Vec3& posB, Contact& contact)
{
    Vec3 delta = posB - posA;
    float distSq = delta.lengthSq();
    float radiusSum = a.radius + b.radius;
    
    if (distSq >= radiusSum * radiusSum) {
        return false; // No collision
    }
    
    float dist = std::sqrt(distSq);
    contact.penetration = radiusSum - dist;
    
    if (dist > 1e-6f) {
        contact.normal = delta / dist;
    } else {
        contact.normal = Vec3(0, 1, 0); // Arbitrary
    }
    
    contact.point = posA + contact.normal * a.radius;
    return true;
}

bool sphereVsBox(const Sphere& sphere, const Vec3& posS, const Quaternion& rotS,
                 const Box& box, const Vec3& posB, const Quaternion& rotB, Contact& contact)
{
    // Transform sphere center to box local space
    Quaternion invRotB = rotB.conjugate();
    Vec3 localSpherePos = invRotB.rotate(posS - posB);
    
    // Closest point on box to sphere
    Vec3 closest;
    closest.x = clamp(localSpherePos.x, -box.halfExtents.x, box.halfExtents.x);
    closest.y = clamp(localSpherePos.y, -box.halfExtents.y, box.halfExtents.y);
    closest.z = clamp(localSpherePos.z, -box.halfExtents.z, box.halfExtents.z);
    
    Vec3 delta = localSpherePos - closest;
    float distSq = delta.lengthSq();
    
    if (distSq >= sphere.radius * sphere.radius) {
        return false;
    }
    
    float dist = std::sqrt(distSq);
    contact.penetration = sphere.radius - dist;
    
    if (dist > 1e-6f) {
        // Transform normal back to world space
        contact.normal = rotB.rotate(delta / dist);
    } else {
        // Sphere center inside box - use closest axis
        Vec3 depths;
        depths.x = box.halfExtents.x - std::abs(localSpherePos.x);
        depths.y = box.halfExtents.y - std::abs(localSpherePos.y);
        depths.z = box.halfExtents.z - std::abs(localSpherePos.z);
        
        if (depths.x < depths.y && depths.x < depths.z) {
            contact.normal = rotB.rotate(Vec3(localSpherePos.x > 0 ? 1.0f : -1.0f, 0, 0));
        } else if (depths.y < depths.z) {
            contact.normal = rotB.rotate(Vec3(0, localSpherePos.y > 0 ? 1.0f : -1.0f, 0));
        } else {
            contact.normal = rotB.rotate(Vec3(0, 0, localSpherePos.z > 0 ? 1.0f : -1.0f));
        }
    }
    
    contact.point = posS - contact.normal * sphere.radius;
    return true;
}

bool boxVsBox(const Box& a, const Vec3& posA, const Quaternion& rotA,
              const Box& b, const Vec3& posB, const Quaternion& rotB, Contact& contact)
{
    // Simplified SAT (Separating Axis Theorem) implementation
    // For production: full SAT with all 15 axes (6 face + 9 edge cross products)
    
    Mat3 rotMatA = rotA.toMat3();
    Mat3 rotMatB = rotB.toMat3();
    
    Vec3 delta = posB - posA;
    
    // Test face normals of A
    for (int i = 0; i < 3; ++i) {
        Vec3 axis(rotMatA.m[i * 3], rotMatA.m[i * 3 + 1], rotMatA.m[i * 3 + 2]);
        
        float projA = a.halfExtents[i];
        float projB = 0;
        for (int j = 0; j < 3; ++j) {
            Vec3 axisB(rotMatB.m[j * 3], rotMatB.m[j * 3 + 1], rotMatB.m[j * 3 + 2]);
            projB += b.halfExtents[j] * std::abs(axis.dot(axisB));
        }
        
        float dist = std::abs(delta.dot(axis));
        if (dist > projA + projB) {
            return false; // Separating axis found
        }
    }
    
    // Test face normals of B
    for (int i = 0; i < 3; ++i) {
        Vec3 axis(rotMatB.m[i * 3], rotMatB.m[i * 3 + 1], rotMatB.m[i * 3 + 2]);
        
        float projB = b.halfExtents[i];
        float projA = 0;
        for (int j = 0; j < 3; ++j) {
            Vec3 axisA(rotMatA.m[j * 3], rotMatA.m[j * 3 + 1], rotMatA.m[j * 3 + 2]);
            projA += a.halfExtents[j] * std::abs(axis.dot(axisA));
        }
        
        float dist = std::abs(delta.dot(axis));
        if (dist > projA + projB) {
            return false;
        }
    }
    
    // Collision detected (simplified contact generation)
    contact.normal = (delta.lengthSq() > 1e-6f) ? delta.normalized() : Vec3(0, 1, 0);
    contact.point = (posA + posB) * 0.5f;
    contact.penetration = 0.1f; // Approximation
    
    return true;
}

} // namespace physx3d
