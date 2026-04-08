#pragma once

// EMBER Integer Exact CSG - Polygon Clipping
// Clips convex polygons against axis-aligned and general planes
// Used during subdivision phase to split polygons at AABB boundaries
//
// Algorithm from paper Section 4.2.1:
//   1. Classify each vertex against splitting plane
//   2. If all on one side, assign polygon to that side
//   3. If straddling, split by inserting splitting plane as new edge

#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>

#include <vector>

namespace ember
{

enum class ClipSide
{
    Left,  // Entirely on negative side of splitting plane
    Right, // Entirely on positive side
    Both   // Straddling - both halves produced
};

struct ClipResult
{
    ConvexPolygon left;  // Polygon on negative side (classify <= 0)
    ConvexPolygon right; // Polygon on positive side (classify >= 0)
    ClipSide side;
};

// Clip a convex polygon against a plane.
// "Left" = negative side of plane (classify <= 0)
// "Right" = positive side (classify >= 0)
inline ClipResult clip_polygon(ConvexPolygon const& poly, plane_t const& split_plane)
{
    int n = poly.vertex_count();
    if (n < 3) return {poly, {}, ClipSide::Left};

    // Step 1: Classify each vertex against the splitting plane
    std::vector<tg::i8> c(n);
    bool has_pos = false, has_neg = false;

    for (int i = 0; i < n; i++)
    {
        auto v = poly.vertex(i);
        c[i] = exact_classify(v, split_plane);
        if (c[i] > 0) has_pos = true;
        if (c[i] < 0) has_neg = true;
    }

    // Step 2: Check for trivial cases
    if (!has_pos)
    {
        // All vertices on negative side or on the plane → entirely left
        return {poly, {}, ClipSide::Left};
    }
    if (!has_neg)
    {
        // All vertices on positive side or on the plane → entirely right
        return {{}, poly, ClipSide::Right};
    }

    // Step 3: Split - walk around the polygon building both halves
    // When transitioning from left→right or right→left, the splitting plane
    // becomes a new edge plane.
    //
    // Convention: vertex i = intersect(support, edges[i], edges[(i+1)%n])
    //             segment v_i→v_{i+1} lies on edges[(i+1)%n]

    plane_t q_inv = split_plane.inverted();

    std::vector<plane_t> left_edges, right_edges;
    left_edges.reserve(n + 2);
    right_edges.reserve(n + 2);

    for (int i = 0; i < n; i++)
    {
        int j = (i + 1) % n;
        plane_t seg_edge = poly.edges[j]; // Edge plane for segment v_i → v_j

        if (c[i] <= 0 && c[j] <= 0)
        {
            // Both on left: segment entirely to left
            left_edges.push_back(seg_edge);
        }
        else if (c[i] <= 0 && c[j] > 0)
        {
            // Transition left → right
            // Split point is intersect(support, seg_edge, split_plane)
            // Left side: include seg_edge (up to split point), then follow split_plane
            // Right side: include seg_edge (from split point onward)
            left_edges.push_back(seg_edge);
            left_edges.push_back(split_plane);
            right_edges.push_back(seg_edge);
        }
        else if (c[i] > 0 && c[j] <= 0)
        {
            // Transition right → left
            // Right side: include seg_edge (up to split point), then follow inverted split_plane
            // Left side: include seg_edge (from split point onward)
            right_edges.push_back(seg_edge);
            right_edges.push_back(q_inv);
            left_edges.push_back(seg_edge);
        }
        else // c[i] > 0 && c[j] > 0
        {
            // Both on right: segment entirely to right
            right_edges.push_back(seg_edge);
        }
    }

    // Construct result polygons
    ClipResult result;
    result.side = ClipSide::Both;

    result.left.support = poly.support;
    result.left.edges = std::move(left_edges);
    result.left.mesh_index = poly.mesh_index;
    result.left.polygon_index = poly.polygon_index;
    result.left.delta_w = poly.delta_w;
    result.left.no_self_intersections = poly.no_self_intersections;
    result.left.no_nested_components = poly.no_nested_components;

    result.right.support = poly.support;
    result.right.edges = std::move(right_edges);
    result.right.mesh_index = poly.mesh_index;
    result.right.polygon_index = poly.polygon_index;
    result.right.delta_w = poly.delta_w;
    result.right.no_self_intersections = poly.no_self_intersections;
    result.right.no_nested_components = poly.no_nested_components;

    return result;
}

// Clip a polygon against an AABB (keep only the part inside the box)
// Returns empty polygon (edges.empty()) if entirely outside
inline ConvexPolygon clip_polygon_to_aabb(ConvexPolygon const& poly, IAABB const& aabb)
{
    ConvexPolygon current = poly;

    // Clip against each of the 6 AABB faces
    // For axis i, min face: normal = -axis, d = min[i]  →  plane: -x[i] + min[i] <= 0  →  x[i] >= min[i]
    // For axis i, max face: normal = +axis, d = -max[i] →  plane: x[i] - max[i] <= 0   →  x[i] <= max[i]
    for (int axis = 0; axis < 3; axis++)
    {
        if (current.edges.empty()) break;

        // Clip against min face: keep x[axis] >= min[axis]
        // Plane: -x[axis] + min[axis] = 0  →  a=-1 (axis=0), b=-1 (axis=1), c=-1 (axis=2)
        // Points with x[axis] >= min[axis] are on the NEGATIVE side of this plane: -x[axis] + min[axis] <= 0
        // So we want to keep the RIGHT side (positive side means x[axis] < min[axis])
        // Actually: let's use plane x[axis] - min[axis] = 0 with normal pointing positive
        // classify > 0 means x[axis] > min[axis]  → we want right side
        {
            plane_t min_plane;
            min_plane.a = (axis == 0) ? normal_scalar_t(1) : normal_scalar_t(0);
            min_plane.b = (axis == 1) ? normal_scalar_t(1) : normal_scalar_t(0);
            min_plane.c = (axis == 2) ? normal_scalar_t(1) : normal_scalar_t(0);
            min_plane.d = plane_d_t(-int64_t((&aabb.min.x)[axis]));

            // We want to keep the positive side (x[axis] >= min[axis])
            // clip_polygon returns left=negative, right=positive
            // So we want the right side
            auto cr = clip_polygon(current, min_plane);
            if (cr.side == ClipSide::Left)
            {
                current.edges.clear();
                break;
            }
            else if (cr.side == ClipSide::Right || cr.side == ClipSide::Both)
            {
                current = cr.side == ClipSide::Right ? current : std::move(cr.right);
            }
        }

        if (current.edges.empty()) break;

        // Clip against max face: keep x[axis] <= max[axis]
        // Plane: x[axis] - max[axis] = 0 with normal pointing positive
        // classify > 0 means x[axis] > max[axis]  → we want left side
        {
            plane_t max_plane;
            max_plane.a = (axis == 0) ? normal_scalar_t(1) : normal_scalar_t(0);
            max_plane.b = (axis == 1) ? normal_scalar_t(1) : normal_scalar_t(0);
            max_plane.c = (axis == 2) ? normal_scalar_t(1) : normal_scalar_t(0);
            max_plane.d = plane_d_t(-int64_t((&aabb.max.x)[axis]));

            auto cr = clip_polygon(current, max_plane);
            if (cr.side == ClipSide::Right)
            {
                current.edges.clear();
                break;
            }
            else if (cr.side == ClipSide::Both)
            {
                current = std::move(cr.left);
            }
            // If Left, current is already correct
        }
    }

    return current;
}

} // namespace ember
