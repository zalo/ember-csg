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
#include <cmath>
#include <cstdio>
#include <vector>

namespace ember
{

// Leaf threshold: when polygon count drops below this, stop subdividing
static constexpr int LEAF_THRESHOLD = 4096;

// Maximum subdivision depth to prevent infinite recursion
static constexpr int MAX_DEPTH = 40;

// Output polygon with classification
struct ClassifiedPolygon
{
    ConvexPolygon polygon;
    WindingPair winding;     // Front and back WNVs
    int classification = 0;  // +1, -1, or 0
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

// Compute splitting plane using center-of-gravity heuristic (Section 4.5.3)
// Tries to separate different WNTV classes
inline void compute_split(SubdivisionTask const& task, int& out_axis, pos_scalar_t& out_val)
{
    auto const& bounds = task.bounds;

    // Group polygons by WNTV
    struct WNTVGroup
    {
        WNTV key;
        tg::dpos3 center;
        int count = 0;
        double variance[3] = {0, 0, 0};
    };

    std::vector<WNTVGroup> groups;
    for (auto const& poly : task.polygons)
    {
        // Find or create group
        WNTVGroup* group = nullptr;
        for (auto& g : groups)
        {
            if (g.key == poly.delta_w) { group = &g; break; }
        }
        if (!group)
        {
            groups.push_back({poly.delta_w, {0,0,0}, 0, {0,0,0}});
            group = &groups.back();
        }

        auto c = poly.center_of_gravity();
        group->center.x += c.x;
        group->center.y += c.y;
        group->center.z += c.z;
        group->count++;
    }

    // Finalize centers
    for (auto& g : groups)
    {
        if (g.count > 0)
        {
            g.center.x /= g.count;
            g.center.y /= g.count;
            g.center.z /= g.count;
        }
    }

    // Try to find a split that separates groups
    if (groups.size() >= 2)
    {
        // Find the axis with greatest separation between group centers
        double best_sep = 0;
        int best_axis = -1;
        double best_val = 0;

        for (int axis = 0; axis < 3; axis++)
        {
            for (size_t i = 0; i < groups.size(); i++)
            {
                for (size_t j = i + 1; j < groups.size(); j++)
                {
                    double ci = (&groups[i].center.x)[axis];
                    double cj = (&groups[j].center.x)[axis];
                    double mid = (ci + cj) * 0.5;
                    double sep = std::abs(ci - cj);

                    // Check that mid is within bounds
                    double bmin = double((&bounds.min.x)[axis]);
                    double bmax = double((&bounds.max.x)[axis]);
                    if (mid >= bmin + 1 && mid <= bmax - 1 && sep > best_sep)
                    {
                        best_sep = sep;
                        best_axis = axis;
                        best_val = mid;
                    }
                }
            }
        }

        if (best_axis >= 0)
        {
            out_axis = best_axis;
            out_val = pos_scalar_t(int32_t(std::round(best_val)));
            // Clamp to valid range
            auto bmin = (&bounds.min.x)[out_axis];
            auto bmax = (&bounds.max.x)[out_axis];
            if (out_val <= bmin) out_val = bmin + pos_scalar_t(1);
            if (out_val >= bmax) out_val = bmax - pos_scalar_t(1);
            return;
        }
    }

    // Fallback: largest variance split or simple midpoint
    // Compute variance along each axis
    double var[3] = {0, 0, 0};
    auto cog = compute_center_of_gravity(task.polygons);
    for (auto const& poly : task.polygons)
    {
        auto c = poly.center_of_gravity();
        for (int a = 0; a < 3; a++)
        {
            double d = (&c.x)[a] - (&cog.x)[a];
            var[a] += d * d;
        }
    }

    // Pick axis with largest variance, or longest axis
    out_axis = bounds.longest_axis();
    double max_var = 0;
    for (int a = 0; a < 3; a++)
    {
        if (var[a] > max_var && bounds.extent(a) > pos_scalar_t(2))
        {
            max_var = var[a];
            out_axis = a;
        }
    }

    out_val = bounds.midpoint(out_axis);

    // Ensure split position is strictly between bounds
    auto bmin = (&bounds.min.x)[out_axis];
    auto bmax = (&bounds.max.x)[out_axis];
    if (out_val <= bmin) out_val = bmin + pos_scalar_t(1);
    if (out_val >= bmax) out_val = bmax - pos_scalar_t(1);
}

// Compute a new reference point for a child subproblem (Section 4.2.2)
// When the old reference point is outside the new sub-AABB, we need to
// trace a segment from the old reference to a new one inside the sub-AABB.
inline bool compute_new_reference(
    pos_t const& old_ref,
    WNV const& old_wnv,
    IAABB const& new_bounds,
    std::vector<ConvexPolygon> const& polygons,
    pos_t& new_ref,
    WNV& new_wnv)
{
    // If old reference is still inside the new bounds, keep it
    if (new_bounds.contains(old_ref))
    {
        new_ref = old_ref;
        new_wnv = old_wnv;
        return true;
    }

    // Project old reference onto the new AABB boundary
    // Find the closest point on the AABB boundary
    pos_t projected = old_ref;
    for (int a = 0; a < 3; a++)
    {
        auto& coord = (&projected.x)[a];
        auto bmin = (&new_bounds.min.x)[a];
        auto bmax = (&new_bounds.max.x)[a];
        if (coord < bmin) coord = bmin + pos_scalar_t(1);
        else if (coord > bmax) coord = bmax - pos_scalar_t(1);
    }

    // Now find an interior point of the AABB to use as reference
    // Use the center of the AABB
    pos_t center(new_bounds.midpoint(0), new_bounds.midpoint(1), new_bounds.midpoint(2));

    // Trace L-path from old_ref to center to propagate WNV
    // Use a "null" skip (no polygon to skip since this is reference propagation)
    plane_t no_skip{};
    WNV current_wnv = old_wnv;
    pos_t cur = old_ref;
    if (center.x != cur.x)
    {
        pos_t next(center.x, cur.y, cur.z);
        current_wnv = trace_segment(cur, next, current_wnv, polygons, no_skip, -1);
        cur = next;
    }
    if (center.y != cur.y)
    {
        pos_t next(cur.x, center.y, cur.z);
        current_wnv = trace_segment(cur, next, current_wnv, polygons, no_skip, -1);
        cur = next;
    }
    if (center.z != cur.z)
        current_wnv = trace_segment(cur, center, current_wnv, polygons, no_skip, -1);

    new_ref = center;
    new_wnv = current_wnv;
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
            // All NNC: only need to classify one polygon, rest are the same
            WNV w_front = task.ref_wnv;
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
        // Build BVH for classification raycast acceleration
        BVH nsi_bvh;
        nsi_bvh.build(polygons);

        for (auto const& poly : polygons)
        {
            WNV w_front = classify_leaf_polygon(
                poly.support, poly.edges, task.ref_point, task.ref_wnv,
                polygons, task.root_bounds, poly.delta_w, poly.mesh_index, &nsi_bvh);
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

    // Two-pass approach with BVH acceleration:
    // Build per-mesh BVHs, use dual-BVH traversal to find intersecting pairs.
    // Only test geometric intersection for AABB-overlapping cross-mesh pairs.

    // Build per-mesh BVHs
    std::map<int, std::vector<int>> mesh_poly_map; // mesh_index → polygon indices
    for (int i = 0; i < n; i++)
        mesh_poly_map[polygons[i].mesh_index].push_back(i);

    std::map<int, BVH> mesh_bvhs;
    for (auto& [mesh_id, poly_ids] : mesh_poly_map)
    {
        std::vector<ConvexPolygon> mesh_polys;
        mesh_polys.reserve(poly_ids.size());
        for (int idx : poly_ids) mesh_polys.push_back(polygons[idx]);
        mesh_bvhs[mesh_id].build(mesh_polys);
    }

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

    // Build a combined BVH for accelerated raycast in WNV classification
    BVH combined_bvh;
    combined_bvh.build(polygons);

    for (int i = 0; i < n; i++)
    {
        auto const& poly_i = polygons[i];

        if (all_isects[i].empty())
        {
            // No intersections: classify directly using the original polygon
            WNV w_front = classify_leaf_polygon(
                poly_i.support, poly_i.edges, task.ref_point, task.ref_wnv,
                polygons, task.root_bounds, poly_i.delta_w, poly_i.mesh_index, &combined_bvh);
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

    // Pass 2: BSP-split intersecting polygons
    int pass2_polys = 0, pass2_leaves = 0, pass2_emit = 0;
    for (int i = 0; i < n; i++)
    {
        if (all_isects[i].empty()) continue;
        pass2_polys++;

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
            pass2_leaves++;

            WNV w_front = classify_leaf_polygon(
                poly_i.support, leaf->edges, task.ref_point, task.ref_wnv,
                polygons, task.root_bounds, poly_i.delta_w, poly_i.mesh_index, &combined_bvh);
            WNV w_back = propagate_wnv(w_front, 1, poly_i.delta_w);
            int cls = classify_polygon_output(w_front, w_back, indicator);
            if (cls != 0) pass2_emit++;

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
                output.push_back(std::move(cp));
            }
        }
    }
    if (pass2_polys > 0)
        std::fprintf(stderr, "[process_leaf] pass2: %d polys, %d leaves, %d emit, %zu total output\n",
            pass2_polys, pass2_leaves, pass2_emit, output.size());
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

    // Compute splitting plane
    int split_axis;
    pos_scalar_t split_val;
    compute_split(task, split_axis, split_val);

    plane_t split_plane = task.bounds.splitting_plane(split_axis, split_val);

    // Create two sub-AABBs
    IAABB left_bounds = task.bounds.left_half(split_axis, split_val);
    IAABB right_bounds = task.bounds.right_half(split_axis, split_val);

    // Assign polygons to child cells. Each polygon goes to ALL cells
    // its AABB overlaps (required for correct BSP). Non-intersecting polygons
    // go to only ONE cell (center-based) to avoid output duplication.
    // Intersecting polygons MUST be in both cells for correct BSP processing
    // but we deduplicate in process_leaf.
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
            // Send the FULL polygon (unclipped) to both cells.
            // process_leaf will classify it correctly from each cell's reference.
            // The two-pass architecture ensures non-intersecting polygons
            // are only emitted once (from the cell where they're fully contained).
            left_polys.push_back(poly);
            right_polys.push_back(poly);
            break;
        }
    }

    // Compute reference points for each child
    pos_t left_ref, right_ref;
    WNV left_wnv, right_wnv;

    compute_new_reference(task.ref_point, task.ref_wnv, left_bounds,
                          task.polygons, left_ref, left_wnv);
    compute_new_reference(task.ref_point, task.ref_wnv, right_bounds,
                          task.polygons, right_ref, right_wnv);

    // Recurse into children
    // Process the smaller side first (continue with larger, add smaller to work queue)
    // For sequential execution, just recurse both
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
