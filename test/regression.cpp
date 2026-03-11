#include <physx3d/core.hpp>
#include <physx3d/spatial.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace physx3d;

static bool finiteVec(const Vec3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static bool check(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "[FAIL] " << msg << "\n";
        return false;
    }
    return true;
}

static bool test_bvh_build_stability() {
    BVH bvh;
    std::vector<AABB> aabbs;
    std::vector<int> indices;

    for (int n = 1; n <= 2048; n += 257) {
        aabbs.clear();
        indices.clear();
        aabbs.reserve(n);
        indices.reserve(n);

        for (int i = 0; i < n; ++i) {
            float x = static_cast<float>((i * 37) % 97) * 0.5f;
            float y = static_cast<float>((i * 53) % 89) * 0.25f;
            float z = static_cast<float>((i * 29) % 83) * 0.4f;
            Vec3 p(x, y, z);
            aabbs.emplace_back(p - Vec3(0.5f, 0.5f, 0.5f), p + Vec3(0.5f, 0.5f, 0.5f));
            indices.push_back(i);
        }

        bvh.build(aabbs, indices);

        std::vector<std::pair<int, int>> pairs;
        bvh.queryOverlaps(pairs);

        int hitIndex = -1;
        bool hit = bvh.raycast(Vec3(0, 0, 0), Vec3(1, 0.1f, 0.2f), 1e6f, hitIndex);
        if (hit && !check(hitIndex >= 0 && hitIndex < n, "BVH raycast returned invalid hit index")) {
            return false;
        }
    }

    return true;
}

static bool test_invalid_mass_and_dt_handling() {
    RigidBody rb;
    rb.setMass(2.0f);
    if (!check(std::abs(rb.getMass() - 2.0f) < 1e-6f, "setMass valid value failed")) return false;

    rb.setMass(-1.0f);
    if (!check(std::abs(rb.getMass() - 2.0f) < 1e-6f, "negative mass was not rejected")) return false;

    rb.setMass(std::numeric_limits<float>::infinity());
    if (!check(std::abs(rb.getMass() - 2.0f) < 1e-6f, "infinite mass was not rejected")) return false;

    rb.setPosition(Vec3(1, 2, 3));
    rb.setVelocity(Vec3(4, 5, 6));
    const Vec3 posBefore = rb.getPosition();
    const Vec3 velBefore = rb.getVelocity();
    rb.integrate(-0.1f);
    if (!check((rb.getPosition() - posBefore).lengthSq() < 1e-12f, "integrate accepted invalid dt (position changed)")) return false;
    if (!check((rb.getVelocity() - velBefore).lengthSq() < 1e-12f, "integrate accepted invalid dt (velocity changed)")) return false;

    PhysicsWorld world;
    auto* body = world.createRigidBody();
    body->setPosition(Vec3(0, 10, 0));
    const Vec3 worldPosBefore = body->getPosition();
    world.step(std::numeric_limits<float>::quiet_NaN());
    if (!check((body->getPosition() - worldPosBefore).lengthSq() < 1e-12f, "world.step accepted NaN dt")) return false;

    return true;
}

static bool test_raycast_hit_miss() {
    PhysicsWorld world(Vec3(0, 0, 0));
    auto* body = world.createRigidBody();
    body->setShape(Sphere(1.0f));
    body->setPosition(Vec3(0, 0, 5));

    RigidBody* hitBody = nullptr;
    Vec3 hitPoint;

    bool hit = world.raycast(Vec3(0, 0, 0), Vec3(0, 0, 1), 100.0f, hitBody, hitPoint);
    if (!check(hit, "raycast expected hit")) return false;
    if (!check(hitBody == body, "raycast returned wrong hit body")) return false;
    if (!check(std::abs(hitPoint.z - 4.0f) < 1e-3f, "raycast returned unexpected hit point")) return false;

    hitBody = nullptr;
    hit = world.raycast(Vec3(0, 0, 0), Vec3(1, 0, 0), 100.0f, hitBody, hitPoint);
    if (!check(!hit, "raycast miss case returned hit")) return false;
    if (!check(hitBody == nullptr, "raycast miss case produced body")) return false;

    return true;
}

static bool stress_many_bodies() {
    PhysicsWorld world(Vec3(0, -9.81f, 0));

    auto* ground = world.createRigidBody();
    ground->setShape(Box(Vec3(200, 1, 200)));
    ground->setPosition(Vec3(0, -1, 0));
    ground->setStatic(true);

    std::vector<RigidBody*> bodies;
    const int bodyCount = 1200;
    const int stepCount = 120;

    bodies.reserve(bodyCount);

    for (int i = 0; i < bodyCount; ++i) {
        auto* b = world.createRigidBody();
        b->setShape(Box(Vec3(0.25f, 0.25f, 0.25f)));
        float x = static_cast<float>((i % 50) - 25) * 0.75f;
        float y = 2.0f + static_cast<float>((i / 50) % 30) * 0.8f;
        float z = static_cast<float>((i / 600) * 3 - 1) * 2.0f;
        b->setPosition(Vec3(x, y, z));
        b->setMass(1.0f);
        b->setRestitution(0.1f);
        bodies.push_back(b);
    }

    const float dt = 1.0f / 120.0f;
    for (int s = 0; s < stepCount; ++s) {
        world.step(dt);
    }

    for (auto* b : bodies) {
        const Vec3 p = b->getPosition();
        const Vec3 v = b->getVelocity();
        if (!check(finiteVec(p), "stress: non-finite position")) return false;
        if (!check(finiteVec(v), "stress: non-finite velocity")) return false;
        if (!check(std::abs(v.x) < 1e4f, "stress: exploding velocity x")) return false;
        if (!check(std::abs(v.y) < 1e4f, "stress: exploding velocity y")) return false;
        if (!check(std::abs(v.z) < 1e4f, "stress: exploding velocity z")) return false;
    }

    return true;
}

int main() {
    if (!test_bvh_build_stability()) return 1;
    if (!test_invalid_mass_and_dt_handling()) return 1;
    if (!test_raycast_hit_miss()) return 1;
    if (!stress_many_bodies()) return 1;

    std::cout << "regression tests passed\n";
    return 0;
}
