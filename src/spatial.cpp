#include "../include/physx3d/spatial.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdint>
#include <cstring>

namespace physx3d {

namespace {

bool isFiniteFloat(float v) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &v, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
}

bool rayIntersectsAABB(const Vec3& origin,
                       const Vec3& direction,
                       const AABB& aabb,
                       float maxDist,
                       float& outT)
{
    float tMin = 0.0f;
    float tMax = maxDist;

    for (int axis = 0; axis < 3; ++axis) {
        float o = origin[axis];
        float d = direction[axis];
        float minA = aabb.min[axis];
        float maxA = aabb.max[axis];

        if (std::abs(d) < 1e-8f) {
            if (o < minA || o > maxA) {
                return false;
            }
            continue;
        }

        float invD = 1.0f / d;
        float t1 = (minA - o) * invD;
        float t2 = (maxA - o) * invD;
        if (t1 > t2) std::swap(t1, t2);

        tMin = std::max(tMin, t1);
        tMax = std::min(tMax, t2);

        if (tMin > tMax) {
            return false;
        }
    }

    outT = tMin;
    return true;
}

} // namespace

void BVH::build(const std::vector<AABB>& aabbs, const std::vector<int>& indices) {
    nodes.clear();
    if (indices.empty()) return;

    buildRecursive(aabbs, indices, 0, static_cast<int>(indices.size()));
}

int BVH::buildRecursive(const std::vector<AABB>& aabbs,
                        const std::vector<int>& indices,
                        int start, int end)
{
    int nodeIdx = static_cast<int>(nodes.size());
    nodes.emplace_back();

    // Compute bounding box for objects referenced by indices[start..end)
    nodes[nodeIdx].aabb = aabbs[indices[start]];
    for (int i = start + 1; i < end; ++i) {
        nodes[nodeIdx].aabb.merge(aabbs[indices[i]]);
    }

    int count = end - start;

    // Leaf node
    if (count == 1) {
        nodes[nodeIdx].objectIndex = indices[start];
        return nodeIdx;
    }

    // Split along longest axis
    Vec3 extent = nodes[nodeIdx].aabb.extents();
    int axis = extent.x > extent.y ? (extent.x > extent.z ? 0 : 2) : (extent.y > extent.z ? 1 : 2);

    // Sort along axis
    std::vector<int> sorted(indices.begin() + start, indices.begin() + end);
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        return aabbs[a].center()[axis] < aabbs[b].center()[axis];
    });

    int mid = count / 2;

    // Build children
    int leftChild = buildRecursive(aabbs, sorted, 0, mid);
    int rightChild = buildRecursive(aabbs, sorted, mid, count);
    nodes[nodeIdx].leftChild = leftChild;
    nodes[nodeIdx].rightChild = rightChild;

    return nodeIdx;
}

void BVH::queryOverlaps(std::vector<std::pair<int, int>>& pairs) const {
    pairs.clear();
    if (nodes.empty()) return;

    queryNode(0, pairs);
}

void BVH::queryNode(int nodeIdx, std::vector<std::pair<int, int>>& pairs) const {
    if (nodeIdx < 0 || nodeIdx >= static_cast<int>(nodes.size())) return;

    const BVHNode& node = nodes[nodeIdx];

    if (node.isLeaf()) return;

    // Check left vs right
    queryNodePair(node.leftChild, node.rightChild, pairs);

    // Recurse
    queryNode(node.leftChild, pairs);
    queryNode(node.rightChild, pairs);
}

void BVH::queryNodePair(int nodeA, int nodeB, std::vector<std::pair<int, int>>& pairs) const {
    if (nodeA < 0 || nodeB < 0) return;

    const BVHNode& a = nodes[nodeA];
    const BVHNode& b = nodes[nodeB];

    if (!a.aabb.overlaps(b.aabb)) return;

    if (a.isLeaf() && b.isLeaf()) {
        int idxA = a.objectIndex;
        int idxB = b.objectIndex;
        if (idxA < idxB) {
            pairs.emplace_back(idxA, idxB);
        } else {
            pairs.emplace_back(idxB, idxA);
        }
        return;
    }

    if (a.isLeaf()) {
        queryNodePair(nodeA, b.leftChild, pairs);
        queryNodePair(nodeA, b.rightChild, pairs);
    } else if (b.isLeaf()) {
        queryNodePair(a.leftChild, nodeB, pairs);
        queryNodePair(a.rightChild, nodeB, pairs);
    } else {
        queryNodePair(a.leftChild, b.leftChild, pairs);
        queryNodePair(a.leftChild, b.rightChild, pairs);
        queryNodePair(a.rightChild, b.leftChild, pairs);
        queryNodePair(a.rightChild, b.rightChild, pairs);
    }
}

bool BVH::raycast(const Vec3& origin, const Vec3& direction, float maxDist, int& hitIndex) const {
    hitIndex = -1;

    if (nodes.empty() || maxDist < 0.0f) {
        return false;
    }

    if (!isFiniteFloat(maxDist) || !isFiniteFloat(origin.x) || !isFiniteFloat(origin.y) ||
        !isFiniteFloat(origin.z) || !isFiniteFloat(direction.x) || !isFiniteFloat(direction.y) ||
        !isFiniteFloat(direction.z))
    {
        return false;
    }

    float dirLenSq = direction.lengthSq();
    if (dirLenSq <= 1e-12f) {
        return false;
    }

    Vec3 dir = direction / std::sqrt(dirLenSq);
    float closestT = maxDist;

    std::vector<int> stack;
    stack.push_back(0);

    while (!stack.empty()) {
        int nodeIdx = stack.back();
        stack.pop_back();

        if (nodeIdx < 0 || nodeIdx >= static_cast<int>(nodes.size())) continue;

        const BVHNode& node = nodes[nodeIdx];
        float t = 0.0f;
        if (!rayIntersectsAABB(origin, dir, node.aabb, closestT, t)) {
            continue;
        }

        if (node.isLeaf()) {
            closestT = t;
            hitIndex = node.objectIndex;
            continue;
        }

        if (node.leftChild >= 0) stack.push_back(node.leftChild);
        if (node.rightChild >= 0) stack.push_back(node.rightChild);
    }

    return hitIndex >= 0;
}

} // namespace physx3d
