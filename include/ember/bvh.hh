#pragma once

// EMBER BVH wrapper using tinybvh (https://github.com/jbikker/tinybvh)
//
// Used for:
//   1. Dual-BVH traversal to find cross-mesh intersecting polygon pairs
//   2. Ray casting for point-in-mesh classification
//
// The BVH operates on double-precision AABBs computed from polygon vertices.
// All actual geometric decisions use EMBER's exact integer arithmetic.

#include <ember/polygon.hh>
#include <functional>
#include <vector>
#include <cmath>

#include "tiny_bvh.h"

namespace ember
{

struct AABB
{
    float min[3] = {1e30f, 1e30f, 1e30f};
    float max[3] = {-1e30f, -1e30f, -1e30f};
    void expand(double x, double y, double z)
    {
        min[0]=std::min(min[0],(float)x); max[0]=std::max(max[0],(float)x);
        min[1]=std::min(min[1],(float)y); max[1]=std::max(max[1],(float)y);
        min[2]=std::min(min[2],(float)z); max[2]=std::max(max[2],(float)z);
    }
    void expand(AABB const& o) {
        for (int a=0;a<3;a++) { min[a]=std::min(min[a],o.min[a]); max[a]=std::max(max[a],o.max[a]); }
    }
    bool intersects(AABB const& o) const {
        return min[0]<=o.max[0]&&max[0]>=o.min[0]&&
               min[1]<=o.max[1]&&max[1]>=o.min[1]&&
               min[2]<=o.max[2]&&max[2]>=o.min[2];
    }
};

struct BVH
{
    tinybvh::BVH tbvh;
    std::vector<AABB> prim_aabbs;
    int prim_count = 0;

    void build(std::vector<ConvexPolygon> const& polys)
    {
        prim_count = (int)polys.size();
        if (prim_count == 0) return;

        // Compute per-polygon AABBs
        prim_aabbs.resize(prim_count);
        for (int i = 0; i < prim_count; i++)
        {
            AABB& box = prim_aabbs[i];
            box = AABB{};
            int nv = polys[i].vertex_count();
            for (int v = 0; v < nv; v++)
            {
                auto dp = polys[i].vertex_dpos(v);
                box.expand(dp.x, dp.y, dp.z);
            }
        }

        // Build using AABB array (6 floats per primitive: minx, miny, minz, maxx, maxy, maxz)
        // tinybvh's BuildAABB expects bvhvec4 pairs: (min.x, min.y, min.z, 0), (max.x, max.y, max.z, 0)
        std::vector<tinybvh::bvhvec4> aabb_data(prim_count * 2);
        for (int i = 0; i < prim_count; i++)
        {
            auto& b = prim_aabbs[i];
            aabb_data[i*2+0] = tinybvh::bvhvec4(b.min[0], b.min[1], b.min[2], 0);
            aabb_data[i*2+1] = tinybvh::bvhvec4(b.max[0], b.max[1], b.max[2], 0);
        }
        tbvh.BuildAABB(aabb_data.data(), prim_count);
    }

    // Query: find all primitives whose AABB overlaps the given box.
    // Uses a simple BVH traversal with AABB-AABB test.
    void query(AABB const& box, std::vector<int>& results) const
    {
        if (prim_count == 0 || tbvh.bvhNode == nullptr) return;
        query_recursive(0, box, results);
    }

    // Dual-BVH: find all (i, j) pairs where AABBs overlap
    void intersect_pairs(BVH const& other,
                          std::function<void(int, int)> const& callback) const
    {
        if (prim_count == 0 || other.prim_count == 0) return;
        if (tbvh.bvhNode == nullptr || other.tbvh.bvhNode == nullptr) return;
        intersect_recursive(0, other, 0, callback);
    }

    // Raycast: find all primitives whose AABB is hit by the ray
    void raycast(double ox, double oy, double oz,
                 double dx, double dy, double dz,
                 std::function<void(int)> const& callback) const
    {
        if (prim_count == 0 || tbvh.bvhNode == nullptr) return;
        float idx = std::abs(dx) > 1e-20f ? 1.0f/(float)dx : 1e20f;
        float idy = std::abs(dy) > 1e-20f ? 1.0f/(float)dy : 1e20f;
        float idz = std::abs(dz) > 1e-20f ? 1.0f/(float)dz : 1e20f;
        raycast_recursive(0, (float)ox,(float)oy,(float)oz, idx,idy,idz, callback);
    }

private:
    static bool aabb_overlap(tinybvh::BVH::BVHNode const& node, AABB const& box)
    {
        return node.aabbMin.x <= box.max[0] && node.aabbMax.x >= box.min[0] &&
               node.aabbMin.y <= box.max[1] && node.aabbMax.y >= box.min[1] &&
               node.aabbMin.z <= box.max[2] && node.aabbMax.z >= box.min[2];
    }

    static bool node_overlap(tinybvh::BVH::BVHNode const& a, tinybvh::BVH::BVHNode const& b)
    {
        return a.aabbMin.x <= b.aabbMax.x && a.aabbMax.x >= b.aabbMin.x &&
               a.aabbMin.y <= b.aabbMax.y && a.aabbMax.y >= b.aabbMin.y &&
               a.aabbMin.z <= b.aabbMax.z && a.aabbMax.z >= b.aabbMin.z;
    }

    void query_recursive(uint32_t nodeIdx, AABB const& box, std::vector<int>& results) const
    {
        auto const& node = tbvh.bvhNode[nodeIdx];
        if (!aabb_overlap(node, box)) return;
        if (node.isLeaf())
        {
            uint32_t first = node.leftFirst, count = node.triCount;
            for (uint32_t i = 0; i < count; i++)
                results.push_back((int)tbvh.primIdx[first + i]);
            return;
        }
        query_recursive(node.leftFirst, box, results);
        query_recursive(node.leftFirst + 1, box, results);
    }

    void intersect_recursive(uint32_t n1, BVH const& other, uint32_t n2,
                              std::function<void(int,int)> const& cb) const
    {
        auto const& a = tbvh.bvhNode[n1];
        auto const& b = other.tbvh.bvhNode[n2];
        if (!node_overlap(a, b)) return;

        bool a_leaf = a.isLeaf(), b_leaf = b.isLeaf();
        if (a_leaf && b_leaf)
        {
            for (uint32_t i = 0; i < a.triCount; i++)
                for (uint32_t j = 0; j < b.triCount; j++)
                {
                    int ai = (int)tbvh.primIdx[a.leftFirst + i];
                    int bj = (int)other.tbvh.primIdx[b.leftFirst + j];
                    if (prim_aabbs[ai].intersects(other.prim_aabbs[bj]))
                        cb(ai, bj);
                }
            return;
        }
        if (a_leaf)
        {
            intersect_recursive(n1, other, b.leftFirst, cb);
            intersect_recursive(n1, other, b.leftFirst + 1, cb);
        }
        else if (b_leaf)
        {
            intersect_recursive(a.leftFirst, other, n2, cb);
            intersect_recursive(a.leftFirst + 1, other, n2, cb);
        }
        else
        {
            intersect_recursive(a.leftFirst, other, n2, cb);
            intersect_recursive(a.leftFirst + 1, other, n2, cb);
        }
    }

    void raycast_recursive(uint32_t nodeIdx,
                            float ox, float oy, float oz,
                            float idx, float idy, float idz,
                            std::function<void(int)> const& cb) const
    {
        auto const& node = tbvh.bvhNode[nodeIdx];
        float t1x=(node.aabbMin.x-ox)*idx, t2x=(node.aabbMax.x-ox)*idx;
        float t1y=(node.aabbMin.y-oy)*idy, t2y=(node.aabbMax.y-oy)*idy;
        float t1z=(node.aabbMin.z-oz)*idz, t2z=(node.aabbMax.z-oz)*idz;
        float tmin=std::max({std::min(t1x,t2x),std::min(t1y,t2y),std::min(t1z,t2z)});
        float tmax=std::min({std::max(t1x,t2x),std::max(t1y,t2y),std::max(t1z,t2z)});
        if (tmax < 0 || tmin > tmax) return;

        if (node.isLeaf())
        {
            for (uint32_t i = 0; i < node.triCount; i++)
                cb((int)tbvh.primIdx[node.leftFirst + i]);
            return;
        }
        raycast_recursive(node.leftFirst, ox,oy,oz, idx,idy,idz, cb);
        raycast_recursive(node.leftFirst + 1, ox,oy,oz, idx,idy,idz, cb);
    }
};

} // namespace ember
