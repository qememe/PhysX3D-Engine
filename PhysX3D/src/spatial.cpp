#include "../include/physx3d/spatial.hpp"
#include <algorithm>

namespace physx3d {

void BVH::build(const std::vector<AABB>& aabbs, const std::vector<int>& indices) {
    nodes.clear();
    if (indices.empty()) return;
    
    buildRecursive(aabbs, indices, 0, indices.size());
}

int BVH::buildRecursive(const std::vector<AABB>& aabbs, 
                        const std::vector<int>& indices,
                        int start, int end)
{
    int nodeIdx = nodes.size();
    nodes.emplace_back();
    BVHNode& node = nodes[nodeIdx];
    
    // Compute bounding box for all objects
    node.aabb = aabbs[indices[start]];
    for (int i = start + 1; i < end; ++i) {
        node.aabb.merge(aabbs[indices[i]]);
    }
    
    int count = end - start;
    
    // Leaf node
    if (count == 1) {
        node.objectIndex = indices[start];
        return nodeIdx;
    }
    
    // Split along longest axis
    Vec3 extent = node.aabb.extents();
    int axis = extent.x > extent.y ? (extent.x > extent.z ? 0 : 2) : (extent.y > extent.z ? 1 : 2);
    
    // Sort along axis
    std::vector<int> sorted = std::vector<int>(indices.begin() + start, indices.begin() + end);
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        return aabbs[a].center()[axis] < aabbs[b].center()[axis];
    });
    
    int mid = count / 2;
    
    // Build children
    node.leftChild = buildRecursive(aabbs, sorted, 0, mid);
    node.rightChild = buildRecursive(aabbs, sorted, mid, count);
    
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
    // Simplified raycast - can be optimized
    hitIndex = -1;
    return false; // TODO: Implement
}

} // namespace physx3d
