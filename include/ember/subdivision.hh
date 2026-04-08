#pragma once

// EMBER Integer Exact CSG - Recursive Adaptive Subdivision
// Per paper Section 4.1-4.2, 4.5:
//   1. Recursively split AABB into sub-AABBs at axis-aligned planes
//   2. Clip polygons against splitting planes
//   3. Propagate reference points via segment tracing
//   4. At leaf nodes (polygon count < threshold), compute local BSPs
//   5. Apply early termination via WNV reachability analysis
//
// Splitting heuristics (Section 4.5.3):
//   - Default: midpoint of longest axis
//   - Advanced: center-of-gravity split to separate WNTV classes
//
// Leaf threshold: 25 polygons (optimal per paper's profiling)

#include <ember/bvh.hh>
#include <ember/clip.hh>
#include <ember/intersect_polygons.hh>
#include <ember/local_bsp.hh>
#include <ember/polygon.hh>
#include <ember/segment_trace.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace ember
{

// Profiling counters (reset before each boolean_operation)
struct SubdivisionProfile {
    double ms_clip = 0;           // Polygon clipping at split planes
    double ms_ref_propagation = 0;// Reference point WNV propagation (includes BVH build)
    double ms_bvh_build = 0;     // BVH construction (subset of ref_propagation + leaf_isect)
    double ms_compute_split = 0;  // Split plane computation
    double ms_leaf_isect = 0;     // Leaf: pairwise intersection detection
    double ms_leaf_bsp = 0;       // Leaf: BSP construction
    double ms_leaf_classify = 0;  // Leaf: polygon classification (segment trace)
    int n_leaves = 0;
    int n_leaf_polys = 0;
    int n_classify_calls = 0;
    int n_subdivide_calls = 0;
};
inline SubdivisionProfile g_profile;
using HiClock = std::chrono::high_resolution_clock;
inline double ms_since(HiClock::time_point t0) {
    return std::chrono::duration<double, std::milli>(HiClock::now() - t0).count();
}

// Leaf threshold: when polygon count drops below this, stop subdividing
static constexpr int LEAF_THRESHOLD = 25;

// Maximum subdivision depth to prevent infinite recursion
static constexpr int MAX_DEPTH = 40;

// Output polygon with classification
struct ClassifiedPolygon
{
    ConvexPolygon polygon;
    WindingPair winding;     // Front and back WNVs
    int classification = 0;  // +1, -1, or 0
    bool is_bsp_fragment = false; // True if this came from BSP splitting
};

// A subproblem in the subdivision tree
struct SubdivisionTask
{
    std::vector<ConvexPolygon> polygons;
    IAABB bounds;
    IAABB root_bounds;       // Top-level AABB (for raycast endpoint, guaranteed outside meshes)
    pos_t ref_point;         // Reference point with known WNV
    WNV ref_wnv;             // WNV at the reference point
    int depth = 0;

};

// Compute the center of gravity of a set of polygons
inline tg::dpos3 compute_center_of_gravity(std::vector<ConvexPolygon> const& polygons)
{
    tg::dpos3 sum = {0, 0, 0};
    int count = 0;
    for (auto const& p : polygons)
    {
        auto c = p.center_of_gravity();
        sum.x += c.x;
        sum.y += c.y;
        sum.z += c.z;
        count++;
    }
    if (count > 0)
    {
        sum.x /= count;
        sum.y /= count;
        sum.z /= count;
    }
    return sum;
}

// Compute splitting plane: midpoint of longest AABB axis.
// This is O(1) per subdivision level. The paper's CoG heuristic (Section 4.5.3)
// gives marginally better splits but costs O(n) per level due to expensive
// vertex_dpos calls in center_of_gravity.
inline void compute_split(SubdivisionTask const& task, int& out_axis, pos_scalar_t& out_val)
{
    auto const& bounds = task.bounds;
    out_axis = bounds.longest_axis();
    out_val = bounds.midpoint(out_axis);

    // Ensure split position is strictly between bounds
    auto bmin = (&bounds.min.x)[out_axis];
    auto bmax = (&bounds.max.x)[out_axis];
    if (out_val <= bmin) out_val = bmin + pos_scalar_t(1);
    if (out_val >= bmax) out_val = bmax - pos_scalar_t(1);
}

// Compute a new reference point for a child subproblem (Section 4.2.2)
//
// Paper approach: project old_ref onto the child AABB boundary. This yields
// a new position reachable by a SINGLE axis-aligned segment. Trace that one
// segment to propagate the WNV. If the segment hits a polygon edge, try
// small perturbations.
inline bool compute_new_reference(
    pos_t const& old_ref,
    WNV const& old_wnv,
    IAABB const& new_bounds,
    std::vector<ConvexPolygon> const& polygons,
    pos_t& new_ref,
    WNV& new_wnv,
    BVH const* bvh = nullptr)
{
    // If old reference is still inside the new bounds, keep it
    if (new_bounds.contains(old_ref))
    {
        new_ref = old_ref;
        new_wnv = old_wnv;
        return true;
    }

    // Project old_ref onto the nearest face of the child AABB.
    // Only ONE coordinate changes → single axis-aligned segment.
    // Find which axis to project along (the one where old_ref is outside).
    int proj_axis = -1;
    for (int a = 0; a < 3; a++)
    {
        auto coord = (&old_ref.x)[a];
        if (coord < (&new_bounds.min.x)[a] || coord > (&new_bounds.max.x)[a])
        {
            proj_axis = a;
            break;
        }
    }

    if (proj_axis < 0)
    {
        // old_ref is inside bounds (shouldn't happen since we checked above)
        new_ref = old_ref;
        new_wnv = old_wnv;
        return true;
    }

    BVH const* bvh_ptr = bvh;

    // Try multiple target positions along the projection axis with small offsets
    // to avoid landing on polygon edges.
    for (int attempt = 0; attempt < 8; attempt++)
    {
        pos_t target = old_ref;

        // Project along the identified axis, with a small offset for retries
        auto coord = (&old_ref.x)[proj_axis];
        auto bmin = (&new_bounds.min.x)[proj_axis];
        auto bmax = (&new_bounds.max.x)[proj_axis];

        if (coord < bmin)
            (&target.x)[proj_axis] = bmin + pos_scalar_t(1 + attempt);
        else
            (&target.x)[proj_axis] = bmax - pos_scalar_t(1 + attempt);

        // Clamp to bounds
        if ((&target.x)[proj_axis] < bmin + pos_scalar_t(1))
            (&target.x)[proj_axis] = bmin + pos_scalar_t(1);
        if ((&target.x)[proj_axis] > bmax - pos_scalar_t(1))
            (&target.x)[proj_axis] = bmax - pos_scalar_t(1);

        // For retries, also perturb the other axes slightly
        if (attempt >= 2)
        {
            for (int a = 0; a < 3; a++)
            {
                if (a == proj_axis) continue;
                auto& c = (&target.x)[a];
                auto lo = (&new_bounds.min.x)[a] + pos_scalar_t(1);
                auto hi = (&new_bounds.max.x)[a] - pos_scalar_t(1);
                c = c + pos_scalar_t(attempt - 1);
                if (c > hi) c = hi;
                if (c < lo) c = lo;
            }
        }

        // Single-axis trace from old_ref to target
        bool valid = true;
        WNV traced = trace_axis_segment(old_ref, target, proj_axis,
                                        old_wnv, polygons, valid);
        if (valid)
        {
            // If we perturbed other axes, trace those too
            if (attempt >= 2)
            {
                pos_t proj_point = old_ref;
                (&proj_point.x)[proj_axis] = (&target.x)[proj_axis];
                // Re-trace: first the projection axis, then the perturbed axes
                traced = old_wnv;
                pos_t cur = old_ref;
                valid = true;
                for (int a = 0; a < 3 && valid; a++)
                {
                    if ((&cur.x)[a] != (&target.x)[a])
                    {
                        pos_t next = cur;
                        (&next.x)[a] = (&target.x)[a];
                        traced = trace_axis_segment(cur, next, a, traced, polygons, valid);
                        cur = next;
                    }
                }
            }

            if (valid)
            {
                new_ref = target;
                new_wnv = traced;
                return true;
            }
        }
    }

    // All attempts failed — fall back to 7-ray at AABB center
    pos_t center(new_bounds.midpoint(0), new_bounds.midpoint(1), new_bounds.midpoint(2));
    new_ref = center;
    new_wnv = point_in_meshes_robust(
        double(center.x), double(center.y), double(center.z),
        polygons, static_cast<int>(old_wnv.size()));
    return true;
}

// Process a leaf task: compute local BSPs, classify polygons, emit results
inline void process_leaf(
    SubdivisionTask const& task,
    IndicatorFn const& indicator,
    std::vector<ClassifiedPolygon>& output)
{
    auto const& polygons = task.polygons;
    int n = static_cast<int>(polygons.size());

    if (n == 0) return;
    g_profile.n_leaves++;
    g_profile.n_leaf_polys += n;

    auto const& class_polys = polygons;

    // Optimization from Section 4.5.1: if all polygons have the same WNTV
    // and NSI flag is set, skip BSP construction
    bool all_same_wntv = true;
    bool all_nsi = true;
    for (int i = 1; i < n; i++)
    {
        if (polygons[i].delta_w != polygons[0].delta_w)
            all_same_wntv = false;
        if (!polygons[i].no_self_intersections)
            all_nsi = false;
    }

    if (all_same_wntv && all_nsi)
    {
        // All polygons have the same WNTV and no self-intersections
        // → no pairwise intersections possible within this group
        // Just classify each polygon directly

        bool all_nnc = true;
        for (auto const& p : polygons)
            if (!p.no_nested_components) all_nnc = false;

        if (all_nnc && n > 0)
        {
            WNV w_front = classify_leaf_polygon(
                polygons[0].support, polygons[0].edges, task.ref_point, task.ref_wnv,
                class_polys, task.bounds, polygons[0].delta_w, polygons[0].mesh_index);
            WNV w_back = propagate_wnv(w_front, 1, polygons[0].delta_w);

            int cls = classify_polygon_output(w_front, w_back, indicator);
            if (cls != 0)
            {
                for (auto const& poly : polygons)
                {
                    ClassifiedPolygon cp;
                    cp.polygon = poly;
                    cp.winding = {w_front, w_back};
                    cp.classification = cls;
                    output.push_back(std::move(cp));
                }
            }
            return;
        }

        // NSI but not NNC: classify each polygon individually
        for (auto const& poly : polygons)
        {
            WNV w_front = classify_leaf_polygon(
                poly.support, poly.edges, task.ref_point, task.ref_wnv,
                class_polys, task.bounds, poly.delta_w, poly.mesh_index);
            WNV w_back = propagate_wnv(w_front, 1, poly.delta_w);

            int cls = classify_polygon_output(w_front, w_back, indicator);
            if (cls != 0)
            {
                ClassifiedPolygon cp;
                cp.polygon = poly;
                cp.winding = {w_front, w_back};
                cp.classification = cls;
                output.push_back(std::move(cp));
            }
        }
        return;
    }

    // Two-pass approach with BVH acceleration
    auto t_isect_start = HiClock::now();

    // Build per-mesh BVHs
    std::map<int, std::vector<int>> mesh_poly_map; // mesh_index → polygon indices
    for (int i = 0; i < n; i++)
        mesh_poly_map[polygons[i].mesh_index].push_back(i);

    auto t_bvh_leaf = HiClock::now();
    std::map<int, BVH> mesh_bvhs;
    for (auto& [mesh_id, poly_ids] : mesh_poly_map)
    {
        std::vector<ConvexPolygon> mesh_polys;
        mesh_polys.reserve(poly_ids.size());
        for (int idx : poly_ids) mesh_polys.push_back(polygons[idx]);
        mesh_bvhs[mesh_id].build(mesh_polys);
    }
    g_profile.ms_bvh_build += ms_since(t_bvh_leaf);

    // Dual-BVH: find cross-mesh AABB-overlapping pairs
    std::vector<std::vector<PairwiseIntersection>> all_isects(n);

    auto mesh_ids = std::vector<int>();
    for (auto& [k, v] : mesh_poly_map) mesh_ids.push_back(k);

    for (size_t mi = 0; mi < mesh_ids.size(); mi++)
    {
        for (size_t mj = mi + 1; mj < mesh_ids.size(); mj++)
        {
            int mid_i = mesh_ids[mi], mid_j = mesh_ids[mj];
            auto& bvh_i = mesh_bvhs[mid_i];
            auto& bvh_j = mesh_bvhs[mid_j];
            auto& ids_i = mesh_poly_map[mid_i];
            auto& ids_j = mesh_poly_map[mid_j];

            bvh_i.intersect_pairs(bvh_j, [&](int local_i, int local_j) {
                int gi = ids_i[local_i]; // global polygon index
                int gj = ids_j[local_j];

                // Geometric intersection test
                auto isect_ij = intersect_polygons(polygons[gi], polygons[gj], gj);
                if (isect_ij.type == PairwiseIntersection::Type::Segment ||
                    isect_ij.type == PairwiseIntersection::Type::Overlap)
                    all_isects[gi].push_back(isect_ij);

                auto isect_ji = intersect_polygons(polygons[gj], polygons[gi], gi);
                if (isect_ji.type == PairwiseIntersection::Type::Segment ||
                    isect_ji.type == PairwiseIntersection::Type::Overlap)
                    all_isects[gj].push_back(isect_ji);
            });
        }
    }

    g_profile.ms_leaf_isect += ms_since(t_isect_start);

    auto t_classify_start = HiClock::now();
    for (int i = 0; i < n; i++)
    {
        auto const& poly_i = polygons[i];

        if (all_isects[i].empty())
        {
            // No intersections: classify directly using the original polygon
            WNV w_front = classify_leaf_polygon(
                poly_i.support, poly_i.edges, task.ref_point, task.ref_wnv,
                class_polys, task.bounds, poly_i.delta_w, poly_i.mesh_index);
            WNV w_back = propagate_wnv(w_front, 1, poly_i.delta_w);
            int cls = classify_polygon_output(w_front, w_back, indicator);

            if (cls != 0)
            {
                ClassifiedPolygon cp;
                cp.polygon = poly_i;
                cp.winding = {w_front, w_back};
                cp.classification = cls;
                output.push_back(std::move(cp));
            }
        }
    }

    g_profile.ms_leaf_classify += ms_since(t_classify_start);

    // Pass 2: BSP-split intersecting polygons
    auto t_bsp_start = HiClock::now();
    for (int i = 0; i < n; i++)
    {
        if (all_isects[i].empty()) continue;

        auto const& poly_i = polygons[i];
        LocalBSP bsp;
        bsp.init(poly_i);

        for (auto& isect : all_isects[i])
        {
            if (isect.type == PairwiseIntersection::Type::Segment)
                bsp.add_segment(isect.segment.v0, isect.segment.v1, isect.segment.split_plane);
            else if (isect.type == PairwiseIntersection::Type::Overlap)
                bsp.add_overlap(polygons[isect.overlap.other_polygon_idx], polygons[isect.overlap.other_polygon_idx].polygon_index);
        }

        std::vector<BSPLeaf*> leaves;
        bsp.collect_leaves(leaves);

        for (auto* leaf : leaves)
        {
            if (!leaf->enabled) continue;
            if (leaf->edges.size() < 3) continue;

            WNV w_front = classify_leaf_polygon(
                poly_i.support, leaf->edges, task.ref_point, task.ref_wnv,
                class_polys, task.bounds, poly_i.delta_w, poly_i.mesh_index);
            WNV w_back = propagate_wnv(w_front, 1, poly_i.delta_w);
            int cls = classify_polygon_output(w_front, w_back, indicator);

            if (cls != 0)
            {
                ClassifiedPolygon cp;
                cp.polygon.support = poly_i.support;
                cp.polygon.edges = leaf->edges;
                cp.polygon.mesh_index = poly_i.mesh_index;
                cp.polygon.polygon_index = poly_i.polygon_index;
                cp.polygon.delta_w = poly_i.delta_w;
                cp.winding = {w_front, w_back};
                cp.classification = cls;
                cp.is_bsp_fragment = true;
                output.push_back(std::move(cp));
            }
        }
    }
    g_profile.ms_leaf_bsp += ms_since(t_bsp_start);
}

// Main recursive subdivision function
inline void subdivide(
    SubdivisionTask task,
    IndicatorFn const& indicator,
    std::vector<ClassifiedPolygon>& output)
{
    // Check for early termination via WNV reachability
    if (!task.polygons.empty())
    {
        std::vector<WNTV> available_wntvs;
        for (auto const& p : task.polygons)
        {
            bool found = false;
            for (auto const& w : available_wntvs)
                if (w == p.delta_w) { found = true; break; }
            if (!found)
                available_wntvs.push_back(p.delta_w);
        }

        if (can_early_terminate(task.ref_wnv, available_wntvs, indicator))
            return;
    }

    // Base case: if polygon count below threshold or max depth reached, process as leaf
    if (static_cast<int>(task.polygons.size()) <= LEAF_THRESHOLD || task.depth >= MAX_DEPTH)
    {
        process_leaf(task, indicator, output);
        return;
    }

    // Also check if AABB is too small to split
    bool can_split = false;
    for (int a = 0; a < 3; a++)
    {
        if (task.bounds.extent(a) > pos_scalar_t(2))
        {
            can_split = true;
            break;
        }
    }
    if (!can_split)
    {
        process_leaf(task, indicator, output);
        return;
    }

    g_profile.n_subdivide_calls++;
    auto t_split = HiClock::now();

    // Compute splitting plane
    int split_axis;
    pos_scalar_t split_val;
    compute_split(task, split_axis, split_val);
    g_profile.ms_compute_split += ms_since(t_split);

    plane_t split_plane = task.bounds.splitting_plane(split_axis, split_val);

    // Create two sub-AABBs
    IAABB left_bounds = task.bounds.left_half(split_axis, split_val);
    IAABB right_bounds = task.bounds.right_half(split_axis, split_val);

    // Assign polygons to child cells. Each polygon goes to ALL cells
    // its AABB overlaps (required for correct BSP). Non-intersecting polygons
    // go to only ONE cell (center-based) to avoid output duplication.
    // Intersecting polygons MUST be in both cells for correct BSP processing
    // but we deduplicate in process_leaf.
    auto t_clip = HiClock::now();
    std::vector<ConvexPolygon> left_polys, right_polys;
    left_polys.reserve(task.polygons.size());
    right_polys.reserve(task.polygons.size());

    for (auto const& poly : task.polygons)
    {
        auto cr = clip_polygon(poly, split_plane);
        switch (cr.side)
        {
        case ClipSide::Left:
            left_polys.push_back(poly);
            break;
        case ClipSide::Right:
            right_polys.push_back(poly);
            break;
        case ClipSide::Both:
            left_polys.push_back(std::move(cr.left));
            right_polys.push_back(std::move(cr.right));
            break;
        }
    }

    g_profile.ms_clip += ms_since(t_clip);

    auto t_ref = HiClock::now();
    pos_t left_ref, right_ref;
    WNV left_wnv, right_wnv;

    compute_new_reference(task.ref_point, task.ref_wnv, left_bounds,
                          task.polygons, left_ref, left_wnv);
    compute_new_reference(task.ref_point, task.ref_wnv, right_bounds,
                          task.polygons, right_ref, right_wnv);

    g_profile.ms_ref_propagation += ms_since(t_ref);

    SubdivisionTask left_task;
    left_task.polygons = std::move(left_polys);
    left_task.bounds = left_bounds;
    left_task.root_bounds = task.root_bounds;
    left_task.ref_point = left_ref;
    left_task.ref_wnv = left_wnv;
    left_task.depth = task.depth + 1;

    SubdivisionTask right_task;
    right_task.polygons = std::move(right_polys);
    right_task.bounds = right_bounds;
    right_task.root_bounds = task.root_bounds;
    right_task.ref_point = right_ref;
    right_task.ref_wnv = right_wnv;
    right_task.depth = task.depth + 1;

    subdivide(std::move(left_task), indicator, output);
    subdivide(std::move(right_task), indicator, output);
}

} // namespace ember
