#include <physx3d/core.hpp>
#include <physx3d/spatial.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

using namespace physx3d;

static bool finiteFloat(float v) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &v, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
}

static bool finiteVec(const Vec3& v) {
    return finiteFloat(v.x) && finiteFloat(v.y) && finiteFloat(v.z);
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

static bool test_stacking_stability() {
    PhysicsWorld world(Vec3(0, -9.81f, 0));

    auto* ground = world.createRigidBody();
    ground->setShape(Box(Vec3(20, 1, 20)));
    ground->setPosition(Vec3(0, -1, 0));
    ground->setStatic(true);

    std::vector<RigidBody*> stack;
    for (int i = 0; i < 6; ++i) {
        auto* b = world.createRigidBody();
        b->setShape(Box(Vec3(0.5f, 0.5f, 0.5f)));
        b->setPosition(Vec3(0, 0.5f + i * 1.02f, 0));
        b->setMass(1.0f);
        b->setRestitution(0.0f);
        stack.push_back(b);
    }

    const float dt = 1.0f / 240.0f;
    for (int s = 0; s < 1200; ++s) {
        world.step(dt);
    }

    for (size_t i = 0; i < stack.size(); ++i) {
        const Vec3 p = stack[i]->getPosition();
        const Vec3 v = stack[i]->getVelocity();
        if (!check(finiteVec(p), "stacking: non-finite position")) return false;
        if (!check(finiteVec(v), "stacking: non-finite velocity")) return false;
        if (!check(p.y > -4.0f, "stacking: body penetrated deeply below ground")) return false;
        if (!check(std::abs(v.y) < 80.0f, "stacking: unstable vertical velocity")) return false;
    }

    return true;
}

struct Snapshot {
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
};

static bool run_many_body_scenario(float dt, int steps, Snapshot& out) {
    PhysicsWorld world(Vec3(0, -9.81f, 0));

    auto* ground = world.createRigidBody();
    ground->setShape(Box(Vec3(100, 1, 100)));
    ground->setPosition(Vec3(0, -1, 0));
    ground->setStatic(true);

    auto* wallL = world.createRigidBody();
    wallL->setShape(Box(Vec3(1, 40, 100)));
    wallL->setPosition(Vec3(-40, 20, 0));
    wallL->setStatic(true);

    auto* wallR = world.createRigidBody();
    wallR->setShape(Box(Vec3(1, 40, 100)));
    wallR->setPosition(Vec3(40, 20, 0));
    wallR->setStatic(true);

    auto* wallF = world.createRigidBody();
    wallF->setShape(Box(Vec3(100, 40, 1)));
    wallF->setPosition(Vec3(0, 20, 40));
    wallF->setStatic(true);

    auto* wallB = world.createRigidBody();
    wallB->setShape(Box(Vec3(100, 40, 1)));
    wallB->setPosition(Vec3(0, 20, -40));
    wallB->setStatic(true);

    std::vector<RigidBody*> bodies;
    const int bodyCount = 300;
    bodies.reserve(bodyCount);

    for (int i = 0; i < bodyCount; ++i) {
        auto* b = world.createRigidBody();

        if (i % 3 == 0) {
            b->setShape(Sphere(0.35f));
        } else {
            b->setShape(Box(Vec3(0.3f, 0.3f, 0.3f)));
        }

        float x = static_cast<float>((i % 45) - 22) * 1.1f;
        float y = 2.0f + static_cast<float>((i / 45) % 20) * 0.9f;
        float z = static_cast<float>(((i / 900) % 3) - 1) * 1.5f;
        b->setPosition(Vec3(x, y, z));
        b->setMass(0.6f + static_cast<float>(i % 7) * 0.12f);
        b->setRestitution((i % 5 == 0) ? 0.2f : 0.0f);
        bodies.push_back(b);
    }

    for (int s = 0; s < steps; ++s) {
        world.step(dt);
    }

    out.positions.clear();
    out.velocities.clear();
    out.positions.reserve(bodyCount);
    out.velocities.reserve(bodyCount);

    for (auto* b : bodies) {
        const Vec3 p = b->getPosition();
        const Vec3 v = b->getVelocity();

        if (!check(finiteVec(p), "many-body: non-finite position")) return false;
        if (!check(finiteVec(v), "many-body: non-finite velocity")) return false;
        if (!check(std::abs(v.x) < 1500.0f && std::abs(v.y) < 1500.0f && std::abs(v.z) < 1500.0f,
                   "many-body: exploding velocity")) return false;
        if (!check(std::abs(p.x) < 5000.0f && std::abs(p.y) < 5000.0f && std::abs(p.z) < 5000.0f,
                   "many-body: body left world bounds")) return false;

        out.positions.push_back(p);
        out.velocities.push_back(v);
    }

    return true;
}

static bool test_many_body_long_run_and_determinism() {
    Snapshot s1;
    Snapshot s2;

    // Thousands of steps long-run validation.
    if (!run_many_body_scenario(1.0f / 240.0f, 2000, s1)) return false;

    // Determinism: same inputs -> same outputs.
    if (!run_many_body_scenario(1.0f / 240.0f, 2000, s2)) return false;

    if (!check(s1.positions.size() == s2.positions.size(), "determinism: size mismatch")) return false;

    for (size_t i = 0; i < s1.positions.size(); ++i) {
        const Vec3 dp = s1.positions[i] - s2.positions[i];
        const Vec3 dv = s1.velocities[i] - s2.velocities[i];
        if (!check(dp.lengthSq() < 1e-12f, "determinism: position mismatch")) return false;
        if (!check(dv.lengthSq() < 1e-12f, "determinism: velocity mismatch")) return false;
    }

    return true;
}

static bool test_large_timestep_behavior() {
    Snapshot s;
    return run_many_body_scenario(1.0f / 30.0f, 400, s);
}

int main() {
    if (!test_bvh_build_stability()) return 1;
    if (!test_invalid_mass_and_dt_handling()) return 1;
    if (!test_raycast_hit_miss()) return 1;
    if (!test_stacking_stability()) return 1;
    if (!test_many_body_long_run_and_determinism()) return 1;
    if (!test_large_timestep_behavior()) return 1;

    std::cout << "regression tests passed\n";
    return 0;
}
