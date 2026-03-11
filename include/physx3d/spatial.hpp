#pragma once

#include "collision.hpp"
#include <vector>
#include <memory>

namespace physx3d {

// ============================================================================
// BVH Node - Binary tree for spatial partitioning
// ============================================================================
struct BVHNode {
    AABB aabb;
    int leftChild = -1;   // Index of left child (-1 if leaf)
    int rightChild = -1;  // Index of right child
    int objectIndex = -1; // Index of object if leaf node
    
    bool isLeaf() const { return objectIndex != -1; }
};

// ============================================================================
// BVH - Bounding Volume Hierarchy for broad-phase collision detection
// ============================================================================
class BVH {
public:
    BVH() = default;
    
    // Build BVH from AABBs and their indices
    void build(const std::vector<AABB>& aabbs, const std::vector<int>& indices);
    
    // Query overlapping pairs
    void queryOverlaps(std::vector<std::pair<int, int>>& pairs) const;
    
    // Raycast query against BVH node/object AABBs (broad-phase approximation).
    bool raycast(const Vec3& origin, const Vec3& direction, float maxDist, int& hitIndex) const;
    
    void clear() { nodes.clear(); }
    
private:
    std::vector<BVHNode> nodes;
    
    int buildRecursive(const std::vector<AABB>& aabbs, 
                       const std::vector<int>& indices,
                       int start, int end);
    
    void queryNode(int nodeIdx, std::vector<std::pair<int, int>>& pairs) const;
    void queryNodePair(int nodeA, int nodeB, std::vector<std::pair<int, int>>& pairs) const;
};

} // namespace physx3d
