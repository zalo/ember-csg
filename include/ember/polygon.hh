#pragma once

// EMBER Integer Exact CSG - Convex Polygon Representation
// Polygons are represented in plane-based form:
//   supporting plane s + edge planes e_1, ..., e_n
// Vertices are computed as intersect(s, e_i, e_{(i+1)%n})
// Interior: classify(x, s) = 0 AND classify(x, e_i) <= 0 for all i

#include <ember/types.hh>
#include <ember/exact_classify.hh>
#include <ember/winding.hh>

#include <cmath>
#include <cstdint>
#include <vector>

namespace ember
{

struct ConvexPolygon
{
    plane_t support;                // Supporting plane (polygon lies on this plane)
    std::vector<plane_t> edges;     // Edge planes (interior is on negative side of each)

    int mesh_index = -1;            // Which input mesh this polygon belongs to
    int polygon_index = -1;         // Original polygon index (for overlap resolution)
    WNTV delta_w;                   // Winding number transition vector

    // Additional flags from paper Section 4.5.1
    bool no_self_intersections = false; // NSI flag
    bool no_nested_components = false;  // NNC flag

    // Cached approximate AABB in INTEGER space for fast-path spatial queries.
    // Set during prepare_input from integer vertex positions, updated
    // conservatively during clip_polygon. Allows skipping expensive
    // exact_classify calls when a polygon is clearly on one side of a split.
    int32_t approx_min[3] = {-MAX_COORD, -MAX_COORD, -MAX_COORD};
    int32_t approx_max[3] = { MAX_COORD,  MAX_COORD,  MAX_COORD};
    bool has_approx_bounds = false;

    int vertex_count() const { return static_cast<int>(edges.size()); }

    // Compute vertex i as intersection of support, edges[i], edges[(i+1)%n]
    point4_t vertex(int i) const
    {
        int n = vertex_count();
        return ipg::intersect(support, edges[i], edges[(i + 1) % n]);
    }

    // Compute approximate 3D position of vertex i.
    // Uses binary search with exact classify to find the integer coordinate,
    // avoiding the buggy i256 division in typed-geometry.
    tg::dpos3 vertex_dpos(int i) const
    {
        auto pt = vertex(i);
        if (!pt.is_valid()) return {0, 0, 0};

        // For each axis, find the integer coordinate via binary search.
        // We know the coordinate is bounded by the input AABB (MAX_COORD).
        // classify(pt, axis_plane_at_k) tells us if pt[axis] > k, < k, or == k.
        auto find_coord = [&](int axis) -> double {
            plane_t test_plane{};
            (&test_plane.a)[axis] = normal_scalar_t(1);

            int32_t lo = -MAX_COORD - 1, hi = MAX_COORD + 1;

            // Binary search for the integer part
            while (lo < hi)
            {
                int32_t mid = lo + (hi - lo) / 2;
                test_plane.d = plane_d_t(-int64_t(mid));
                auto c = exact_classify(pt, test_plane);
                if (c > 0)
                    lo = mid + 1; // pt[axis] > mid
                else
                    hi = mid; // pt[axis] <= mid
            }
            // lo == hi == floor(pt[axis]/pt.w) + 1 or the exact value

            // Refine: check classify at lo-1 and lo
            test_plane.d = plane_d_t(-int64_t(lo));
            auto c_lo = exact_classify(pt, test_plane);
            if (c_lo == 0) return double(lo); // exact integer

            // The true value is between lo-1 and lo
            // Use linear interpolation for the fractional part:
            // At lo-1: classify > 0 (pt > lo-1)
            // At lo: classify < 0 (pt < lo)
            // True value = lo - 1 + fraction
            // The fraction can be estimated using double arithmetic on the plane values
            test_plane.d = plane_d_t(-int64_t(lo - 1));
            return double(lo - 1) + 0.5; // approximate to midpoint
        };

        return {find_coord(0), find_coord(1), find_coord(2)};
    }

    // Compute approximate center of gravity (for splitting heuristics)
    tg::dpos3 center_of_gravity() const
    {
        tg::dpos3 sum = {0, 0, 0};
        int n = vertex_count();
        for (int i = 0; i < n; i++)
        {
            auto p = vertex_dpos(i);
            sum.x += p.x;
            sum.y += p.y;
            sum.z += p.z;
        }
        if (n > 0)
        {
            sum.x /= n;
            sum.y /= n;
            sum.z /= n;
        }
        return sum;
    }

    // Invert the polygon (flip normal direction, reverse edge winding)
    ConvexPolygon inverted() const
    {
        ConvexPolygon result;
        result.support = support.inverted();
        result.mesh_index = mesh_index;
        result.polygon_index = polygon_index;
        result.delta_w = delta_w;

        // Reverse and invert edge planes
        int n = vertex_count();
        result.edges.resize(n);
        for (int i = 0; i < n; i++)
            result.edges[i] = edges[n - 1 - i].inverted();

        return result;
    }

    // Check if a point4 is inside the polygon (on supporting plane and inside all edges)
    bool contains_point(point4_t const& pt) const
    {
        if (exact_classify(pt, support) != 0)
            return false;
        for (auto const& e : edges)
            if (exact_classify(pt, e) > 0)
                return false;
        return true;
    }

    // Check if a point4 is strictly inside the polygon (not on any edge)
    bool contains_point_strictly(point4_t const& pt) const
    {
        if (exact_classify(pt, support) != 0)
            return false;
        for (auto const& e : edges)
            if (exact_classify(pt, e) >= 0)
                return false;
        return true;
    }
};

// Create a triangle polygon from three integer positions
// The normal is computed as cross(p1-p0, p2-p0) with GCD reduction
inline ConvexPolygon make_triangle(pos_t p0, pos_t p1, pos_t p2, int mesh_idx, int poly_idx)
{
    ConvexPolygon poly;
    poly.support = plane_t::from_points(p0, p1, p2);
    poly.mesh_index = mesh_idx;
    poly.polygon_index = poly_idx;

    // Edge planes: each edge plane contains the edge and has its normal
    // pointing inward (so classify(interior, edge) < 0)
    //
    // For triangle (p0, p1, p2) with supporting plane s:
    //   edge 0: plane through p0, p1 with normal = cross(s.normal, p1-p0)
    //   edge 1: plane through p1, p2 with normal = cross(s.normal, p2-p1)
    //   edge 2: plane through p2, p0 with normal = cross(s.normal, p0-p2)
    //
    // But we need the normal to point INWARD (classify(interior, edge) <= 0).
    // The inward normal for edge (pa, pb) is: cross(pb - pa, s.normal)
    // This is because the polygon interior is on the left when walking from pa to pb
    // with s.normal pointing "up".

    // Edge planes via from_points: pass the two edge vertices plus a third point
    // offset along the support normal. from_points guarantees 55-bit normals
    // because all input positions are 26-bit.
    auto sn = poly.support.normal();
    int dom = 0;
    {
        auto abs_a = tg::abs(int64_t(sn.x));
        auto abs_b = tg::abs(int64_t(sn.y));
        auto abs_c = tg::abs(int64_t(sn.z));
        if (abs_b > abs_a && abs_b > abs_c) dom = 1;
        else if (abs_c > abs_a && abs_c > abs_b) dom = 2;
    }

    struct EdgeInfo { pos_t a, b, opp; };
    EdgeInfo ev[3] = {{p0, p1, p2}, {p1, p2, p0}, {p2, p0, p1}};

    poly.edges.resize(3);
    for (int i = 0; i < 3; i++)
    {
        // Create a third point off the triangle plane along the support normal
        pos_t off = ev[i].a;
        (&off.x)[dom] = pos_scalar_t(int32_t((&off.x)[dom]) + (int64_t((&sn.x)[dom]) > 0 ? 1 : -1));

        auto edge_plane = plane_t::from_points(ev[i].a, ev[i].b, off);

        // Ensure the opposite vertex is on the negative side (interior convention)
        if (ipg::classify<geometry_t>(ev[i].opp, edge_plane) > 0)
            edge_plane = edge_plane.inverted();

        poly.edges[i] = edge_plane;
    }

    return poly;
}

// Create a quad polygon from four coplanar integer positions (in CCW order)
inline ConvexPolygon make_quad(pos_t p0, pos_t p1, pos_t p2, pos_t p3, int mesh_idx, int poly_idx)
{
    ConvexPolygon poly;
    poly.support = plane_t::from_points(p0, p1, p2);
    poly.mesh_index = mesh_idx;
    poly.polygon_index = poly_idx;

    auto sn = poly.support.normal();
    int dom = 0;
    {
        auto abs_a = tg::abs(int64_t(sn.x));
        auto abs_b = tg::abs(int64_t(sn.y));
        auto abs_c = tg::abs(int64_t(sn.z));
        if (abs_b > abs_a && abs_b > abs_c) dom = 1;
        else if (abs_c > abs_a && abs_c > abs_b) dom = 2;
    }

    // For quad, opposite vertex for edge i is the vertex across from it
    pos_t all_pts[4] = {p0, p1, p2, p3};
    // edge i goes from pt[i] to pt[(i+1)%4], opposite is pt[(i+2)%4]

    poly.edges.resize(4);
    for (int i = 0; i < 4; i++)
    {
        pos_t ea = all_pts[i];
        pos_t eb = all_pts[(i + 1) % 4];
        pos_t opp = all_pts[(i + 2) % 4];

        pos_t off = ea;
        (&off.x)[dom] = pos_scalar_t(int32_t((&off.x)[dom]) + (int64_t((&sn.x)[dom]) > 0 ? 1 : -1));

        auto edge_plane = plane_t::from_points(ea, eb, off);
        if (ipg::classify<geometry_t>(opp, edge_plane) > 0)
            edge_plane = edge_plane.inverted();
        poly.edges[i] = edge_plane;
    }

    return poly;
}

} // namespace ember
