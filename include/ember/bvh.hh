#pragma once

// Simple AABB BVH for accelerating polygon intersection queries.
// Based on the construction approach from three-mesh-bvh:
//   - Recursive top-down build with SAH binning (32 bins)
//   - Flat array storage (left child = node+1, right child = offset)
//   - Supports: AABB query, dual-BVH traversal, ray casting

#include <ember/polygon.hh>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace ember
{

struct AABB
{
    double min[3] = {1e30, 1e30, 1e30};
    double max[3] = {-1e30, -1e30, -1e30};

    void expand(double x, double y, double z)
    {
        min[0] = std::min(min[0], x); max[0] = std::max(max[0], x);
        min[1] = std::min(min[1], y); max[1] = std::max(max[1], y);
        min[2] = std::min(min[2], z); max[2] = std::max(max[2], z);
    }
    void expand(AABB const& o)
    {
        for (int a = 0; a < 3; a++) { min[a] = std::min(min[a], o.min[a]); max[a] = std::max(max[a], o.max[a]); }
    }
    bool intersects(AABB const& o) const
    {
        return min[0] <= o.max[0] && max[0] >= o.min[0] &&
               min[1] <= o.max[1] && max[1] >= o.min[1] &&
               min[2] <= o.max[2] && max[2] >= o.min[2];
    }
    double surface_area() const
    {
        double dx = max[0]-min[0], dy = max[1]-min[1], dz = max[2]-min[2];
        return 2.0 * (dx*dy + dy*dz + dz*dx);
    }
    double center(int axis) const { return (min[axis] + max[axis]) * 0.5; }
};

struct BVHNode
{
    AABB bounds;
    int left = -1;       // left child index (-1 = leaf)
    int right = -1;      // right child index
    int prim_start = 0;  // first primitive (for leaves)
    int prim_count = 0;  // number of primitives (>0 = leaf)
};

struct BVH
{
    std::vector<BVHNode> nodes;
    std::vector<int> prim_indices;  // reordered primitive indices
    std::vector<AABB> prim_aabbs;  // per-primitive AABBs
    std::vector<double> prim_centers; // prim_centers[i*3+axis] = centroid[axis]

    // Build BVH over polygon AABBs
    void build(std::vector<ConvexPolygon> const& polys)
    {
        int n = (int)polys.size();
        if (n == 0) return;

        prim_aabbs.resize(n);
        prim_centers.resize(n * 3);
        prim_indices.resize(n);

        for (int i = 0; i < n; i++)
        {
            prim_indices[i] = i;
            AABB box;
            int nv = polys[i].vertex_count();
            for (int v = 0; v < nv; v++)
            {
                auto dp = polys[i].vertex_dpos(v);
                box.expand(dp.x, dp.y, dp.z);
            }
            prim_aabbs[i] = box;
            prim_centers[i*3+0] = box.center(0);
            prim_centers[i*3+1] = box.center(1);
            prim_centers[i*3+2] = box.center(2);
        }

        nodes.reserve(2 * n);
        build_recursive(0, n, 0);
    }

    // Query: find all primitives whose AABB overlaps the given box
    void query(AABB const& box, std::vector<int>& results) const
    {
        if (nodes.empty()) return;
        query_recursive(0, box, results);
    }

    // Dual-BVH: find all pairs (i, j) where prim_aabbs[i] from this BVH
    // overlaps prim_aabbs[j] from the other BVH.
    // Calls callback(prim_index_this, prim_index_other) for each pair.
    void intersect_pairs(BVH const& other,
                          std::function<void(int, int)> const& callback) const
    {
        if (nodes.empty() || other.nodes.empty()) return;
        intersect_recursive(0, other, 0, callback);
    }

    // Ray cast: find all primitives whose AABB is hit by the ray
    void raycast(double ox, double oy, double oz,
                 double dx, double dy, double dz,
                 std::function<void(int)> const& callback) const
    {
        if (nodes.empty()) return;
        double inv_dx = (std::abs(dx) > 1e-20) ? 1.0/dx : 1e20;
        double inv_dy = (std::abs(dy) > 1e-20) ? 1.0/dy : 1e20;
        double inv_dz = (std::abs(dz) > 1e-20) ? 1.0/dz : 1e20;
        raycast_recursive(0, ox, oy, oz, inv_dx, inv_dy, inv_dz, callback);
    }

private:
    static constexpr int MAX_LEAF = 8;
    static constexpr int NUM_BINS = 32;

    int build_recursive(int start, int end, int depth)
    {
        int count = end - start;
        int node_idx = (int)nodes.size();
        nodes.push_back({});
        auto& node = nodes[node_idx];

        // Compute bounds
        for (int i = start; i < end; i++)
            node.bounds.expand(prim_aabbs[prim_indices[i]]);

        if (count <= MAX_LEAF || depth >= 40)
        {
            node.prim_start = start;
            node.prim_count = count;
            return node_idx;
        }

        // SAH binning: find best split
        int best_axis = -1;
        double best_cost = 1e30;
        int best_split = -1;

        double parent_sa = node.bounds.surface_area();
        if (parent_sa < 1e-20)
        {
            node.prim_start = start;
            node.prim_count = count;
            return node_idx;
        }

        for (int axis = 0; axis < 3; axis++)
        {
            double axis_min = 1e30, axis_max = -1e30;
            for (int i = start; i < end; i++)
            {
                double c = prim_centers[prim_indices[i]*3 + axis];
                axis_min = std::min(axis_min, c);
                axis_max = std::max(axis_max, c);
            }
            if (axis_max - axis_min < 1e-10) continue;

            double bin_width = (axis_max - axis_min) / NUM_BINS;
            int bin_count[NUM_BINS] = {};
            AABB bin_bounds[NUM_BINS];

            for (int i = start; i < end; i++)
            {
                double c = prim_centers[prim_indices[i]*3 + axis];
                int bin = std::clamp(int((c - axis_min) / bin_width), 0, NUM_BINS - 1);
                bin_count[bin]++;
                bin_bounds[bin].expand(prim_aabbs[prim_indices[i]]);
            }

            // Sweep from left
            AABB left_box;
            int left_count = 0;
            double left_sa[NUM_BINS], left_cnt[NUM_BINS];
            for (int i = 0; i < NUM_BINS; i++)
            {
                left_count += bin_count[i];
                left_box.expand(bin_bounds[i]);
                left_sa[i] = left_box.surface_area();
                left_cnt[i] = left_count;
            }

            // Sweep from right
            AABB right_box;
            int right_count = 0;
            for (int i = NUM_BINS - 1; i > 0; i--)
            {
                right_count += bin_count[i];
                right_box.expand(bin_bounds[i]);
                double cost = 1.0 + 1.25 * (left_sa[i-1] * left_cnt[i-1] + right_box.surface_area() * right_count) / parent_sa;
                if (cost < best_cost && left_cnt[i-1] > 0 && right_count > 0)
                {
                    best_cost = cost;
                    best_axis = axis;
                    best_split = i;
                }
            }
        }

        if (best_axis < 0)
        {
            node.prim_start = start;
            node.prim_count = count;
            return node_idx;
        }

        // Partition primitives
        double axis_min = 1e30, axis_max = -1e30;
        for (int i = start; i < end; i++)
        {
            double c = prim_centers[prim_indices[i]*3 + best_axis];
            axis_min = std::min(axis_min, c);
            axis_max = std::max(axis_max, c);
        }
        double bin_width = (axis_max - axis_min) / NUM_BINS;
        double split_pos = axis_min + best_split * bin_width;

        int mid = start;
        for (int i = start; i < end; i++)
        {
            if (prim_centers[prim_indices[i]*3 + best_axis] < split_pos)
                std::swap(prim_indices[i], prim_indices[mid++]);
        }
        if (mid == start || mid == end) mid = (start + end) / 2;

        node.left = build_recursive(start, mid, depth + 1);
        node.right = build_recursive(mid, end, depth + 1);
        return node_idx;
    }

    void query_recursive(int node_idx, AABB const& box, std::vector<int>& results) const
    {
        auto const& node = nodes[node_idx];
        if (!node.bounds.intersects(box)) return;
        if (node.prim_count > 0)
        {
            for (int i = 0; i < node.prim_count; i++)
                results.push_back(prim_indices[node.prim_start + i]);
            return;
        }
        query_recursive(node.left, box, results);
        query_recursive(node.right, box, results);
    }

    void intersect_recursive(int n1, BVH const& other, int n2,
                              std::function<void(int,int)> const& cb) const
    {
        auto const& a = nodes[n1];
        auto const& b = other.nodes[n2];
        if (!a.bounds.intersects(b.bounds)) return;

        if (a.prim_count > 0 && b.prim_count > 0)
        {
            for (int i = 0; i < a.prim_count; i++)
                for (int j = 0; j < b.prim_count; j++)
                    if (prim_aabbs[prim_indices[a.prim_start+i]].intersects(
                        other.prim_aabbs[other.prim_indices[b.prim_start+j]]))
                        cb(prim_indices[a.prim_start+i], other.prim_indices[b.prim_start+j]);
            return;
        }

        if (a.prim_count > 0)
        {
            intersect_recursive(n1, other, b.left, cb);
            intersect_recursive(n1, other, b.right, cb);
        }
        else if (b.prim_count > 0)
        {
            intersect_recursive(a.left, other, n2, cb);
            intersect_recursive(a.right, other, n2, cb);
        }
        else
        {
            intersect_recursive(a.left, other, n2, cb);
            intersect_recursive(a.right, other, n2, cb);
        }
    }

    void raycast_recursive(int node_idx,
                            double ox, double oy, double oz,
                            double idx, double idy, double idz,
                            std::function<void(int)> const& cb) const
    {
        auto const& node = nodes[node_idx];
        // Ray-AABB intersection (slab method)
        double t1x = (node.bounds.min[0] - ox) * idx, t2x = (node.bounds.max[0] - ox) * idx;
        double t1y = (node.bounds.min[1] - oy) * idy, t2y = (node.bounds.max[1] - oy) * idy;
        double t1z = (node.bounds.min[2] - oz) * idz, t2z = (node.bounds.max[2] - oz) * idz;
        double tmin = std::max({std::min(t1x,t2x), std::min(t1y,t2y), std::min(t1z,t2z)});
        double tmax = std::min({std::max(t1x,t2x), std::max(t1y,t2y), std::max(t1z,t2z)});
        if (tmax < 0 || tmin > tmax) return;

        if (node.prim_count > 0)
        {
            for (int i = 0; i < node.prim_count; i++)
                cb(prim_indices[node.prim_start + i]);
            return;
        }
        raycast_recursive(node.left, ox, oy, oz, idx, idy, idz, cb);
        raycast_recursive(node.right, ox, oy, oz, idx, idy, idz, cb);
    }
};

} // namespace ember
