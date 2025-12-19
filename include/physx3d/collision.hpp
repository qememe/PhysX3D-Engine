#pragma once

#include "math.hpp"
#include <memory>
#include <vector>
#include <variant>

namespace physx3d {

// ============================================================================
// Shape Types
// ============================================================================

struct Sphere {
    float radius;
    explicit Sphere(float r = 1.0f) : radius(r) {}
};

struct Box {
    Vec3 halfExtents; // Half-size along each axis
    explicit Box(const Vec3& he) : halfExtents(he) {}
};

struct Capsule {
    float radius;
    float halfHeight; // Half-height of cylindrical part
    explicit Capsule(float r = 0.5f, float hh = 1.0f) : radius(r), halfHeight(hh) {}
};

// Shape variant
using Shape = std::variant<Sphere, Box, Capsule>;

// ============================================================================
// AABB - Axis-Aligned Bounding Box
// ============================================================================
struct AABB {
    Vec3 min, max;
    
    AABB() : min(0, 0, 0), max(0, 0, 0) {}
    AABB(const Vec3& min, const Vec3& max) : min(min), max(max) {}
    
    Vec3 center() const { return (min + max) * 0.5f; }
    Vec3 extents() const { return (max - min) * 0.5f; }
    
    bool overlaps(const AABB& other) const {
        return (min.x <= other.max.x && max.x >= other.min.x) &&
               (min.y <= other.max.y && max.y >= other.min.y) &&
               (min.z <= other.max.z && max.z >= other.min.z);
    }
    
    void expand(const Vec3& point) {
        min.x = std::min(min.x, point.x);
        min.y = std::min(min.y, point.y);
        min.z = std::min(min.z, point.z);
        max.x = std::max(max.x, point.x);
        max.y = std::max(max.y, point.y);
        max.z = std::max(max.z, point.z);
    }
    
    void merge(const AABB& other) {
        expand(other.min);
        expand(other.max);
    }
    
    float surfaceArea() const {
        Vec3 d = max - min;
        return 2.0f * (d.x * d.y + d.y * d.z + d.z * d.x);
    }
};

// ============================================================================
// Contact - collision contact point
// ============================================================================
struct Contact {
    Vec3 point;          // World-space contact point
    Vec3 normal;         // Contact normal (from A to B)
    float penetration;   // Penetration depth
    int bodyA, bodyB;    // Body indices
    
    Contact() : penetration(0), bodyA(-1), bodyB(-1) {}
};

// ============================================================================
// Collision Detection Functions
// ============================================================================

// Compute AABB for shape in world space
AABB computeAABB(const Shape& shape, const Vec3& position, const Quaternion& rotation);

// Narrow-phase collision detection
bool detectCollision(
    const Shape& shapeA, const Vec3& posA, const Quaternion& rotA,
    const Shape& shapeB, const Vec3& posB, const Quaternion& rotB,
    Contact& contact
);

// Specific shape-shape tests
bool sphereVsSphere(const Sphere& a, const Vec3& posA,
                    const Sphere& b, const Vec3& posB, Contact& contact);

bool sphereVsBox(const Sphere& sphere, const Vec3& posS, const Quaternion& rotS,
                 const Box& box, const Vec3& posB, const Quaternion& rotB, Contact& contact);

bool boxVsBox(const Box& a, const Vec3& posA, const Quaternion& rotA,
              const Box& b, const Vec3& posB, const Quaternion& rotB, Contact& contact);

} // namespace physx3d
