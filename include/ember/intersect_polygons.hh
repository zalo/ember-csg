#pragma once

// EMBER Integer Exact CSG - Polygon-Polygon Intersection
// Computes pairwise intersection segments between convex polygons
// Cases from paper Section 4.3:
//   C1: No intersection
//   C2: Single point intersection (ignored)
//   C3: Non-degenerate segment
//   C4: Non-empty overlap (coplanar)

#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>

#include <vector>

namespace ember
{

// Intersection segment between two polygons
// Defined by two endpoints (as point4 homogeneous coords) and a splitting plane
struct IntersectionSegment
{
    point4_t v0, v1;       // Segment endpoints
    plane_t split_plane;   // The supporting plane of the "other" polygon
                           // (together with the host polygon's support, defines the intersection line)
    int other_polygon_idx; // Index of the other polygon in the local list
    bool is_valid = false;
};

// Overlap region: two coplanar polygons overlap
struct OverlapInfo
{
    int other_polygon_idx;
    std::vector<plane_t> other_edges; // Edge planes of the other polygon
    plane_t other_support;            // Supporting plane of the other polygon
};

// Result of intersecting polygon t with polygon t'
struct PairwiseIntersection
{
    enum class Type
    {
        None,    // C1: No intersection
        Point,   // C2: Single point (ignored)
        Segment, // C3: Non-degenerate segment
        Overlap  // C4: Coplanar overlap
    };

    Type type = Type::None;
    IntersectionSegment segment; // Valid if type == Segment
    OverlapInfo overlap;         // Valid if type == Overlap
};

// Compute the intersection of a line (defined by two planes) with a convex polygon.
// Returns a segment [v0, v1] if the line intersects the polygon interior,
// or an invalid segment if no intersection.
//
// The line is the intersection of planes p0 and p1.
// The polygon is defined by supporting plane s and edge planes e[0..n-1].
//
// We clip the infinite line against each edge plane to find the bounded segment.
inline IntersectionSegment clip_line_to_polygon(plane_t const& line_plane0,
                                                 plane_t const& line_plane1,
                                                 ConvexPolygon const& poly)
{
    IntersectionSegment result;
    result.is_valid = false;

    // The line is intersect(line_plane0, line_plane1)
    auto line = ipg::intersect<geometry_t>(line_plane0, line_plane1);

    // Check if line is parallel to supporting plane
    if (ipg::are_parallel<geometry_t>(poly.support, line))
        return result;

    // Find intersection of line with supporting plane → point on the polygon's plane
    // But we need a segment, not a point. The intersection of a line with a plane
    // is only a point if the line is not on the plane.
    // For segment computation, we clip the line against each edge plane.

    // We'll track the valid interval along the line using two "boundary" planes.
    // Initially unbounded. Each edge plane clips the interval.

    int n = poly.vertex_count();

    // For each edge plane, classify the line direction against it.
    // The line parameterized as: any_point(line) + t * direction(line)
    // We find t values where the line crosses each edge plane.
    // The valid interval is the intersection of all half-spaces.

    // Instead of parametric, use the plane-based approach:
    // A point on the line is inside the polygon if:
    //   classify(point, support) = 0  AND  classify(point, edge_i) <= 0 for all i
    //
    // The line intersect(line_plane0, line_plane1) lies on both line_plane0 and line_plane1.
    // We need to find the sub-segment also on the support plane → that's just one point
    // (unless the line lies on the support plane).
    //
    // Wait, that's only true if the line is NOT on the support plane.
    // If the line is on the support plane, we get a line segment (the intersection
    // of the line with the polygon).

    // Check if line lies on supporting plane
    // Line lies on support iff support is parallel to line
    // Already checked above - if parallel, returned.
    // So the line intersects the support plane at exactly one point.
    // This means C3 (segment) requires the line to be ON the support plane.

    // Actually, for polygon-polygon intersection:
    // The intersection of two non-parallel polygons is a line defined by
    // the intersection of their supporting planes. This line lies on BOTH
    // supporting planes. We then clip this line against both polygons' edge
    // planes to get the actual intersection segment.

    // So the correct approach for polygon-polygon intersection:
    // 1. Check if supporting planes are parallel → if so, check overlap (C4)
    // 2. Compute intersection line = intersect(support_t, support_t')
    // 3. Clip line against edge planes of polygon t → segment s_t
    // 4. Clip line against edge planes of polygon t' → segment s_t'
    // 5. Intersect s_t and s_t' → final segment

    // This function assumes the line already lies on the polygon's supporting plane.
    // It clips the line against the polygon's edge planes.

    // Track the valid interval using "left" and "right" boundary planes
    // A point x on the line is valid if classify(x, left_bound) >= 0
    //                              AND classify(x, right_bound) >= 0
    // (or equivalently, we track the two "most restrictive" edge planes)

    bool has_left = false, has_right = false;
    plane_t left_bound, right_bound;

    for (int i = 0; i < n; i++)
    {
        auto const& edge = poly.edges[i];

        // Check if line is parallel to this edge
        if (ipg::are_parallel<geometry_t>(edge, line))
        {
            // Line is parallel to edge. Check which side.
            // Get any point on the line
            auto p = ipg::any_point<geometry_t>(line);
            auto c = exact_classify(p, edge);
            if (c > 0)
                return result; // Line is entirely outside this edge → no intersection
            // c <= 0: line is inside or on this edge → edge doesn't constrain
            continue;
        }

        // Line crosses this edge plane. Determine the "direction" of the constraint.
        // We need to know: as we move along the line, do we enter or exit
        // the negative half-space of this edge?
        //
        // Use the dot product of the line direction with the edge normal to determine
        // the sign. If dot > 0, the line moves from negative to positive (exit),
        // so this edge provides a "left" (entry) boundary.
        // If dot < 0, line moves from positive to negative (entry from the other end),
        // so this edge provides a "right" boundary.

        auto line_dir = line.direction();
        auto edge_n = edge.normal();

        // Compute dot(line_dir, edge_normal)
        // line_dir components are nn_t (2*bits_normal+1 bits)
        // edge_n components are normal_scalar_t (bits_normal bits)
        // dot product fits in 3*bits_normal+3 bits
        static constexpr int dot_bits = 3 * geometry_t::bits_normal + 3;
        auto dot = ipg::mul<dot_bits>(line_dir.x, edge_n.x) +
                   ipg::mul<dot_bits>(line_dir.y, edge_n.y) +
                   ipg::mul<dot_bits>(line_dir.z, edge_n.z);

        if (tg::sign(dot) > 0)
        {
            // Line enters negative half-space from positive: left boundary
            if (!has_left)
            {
                left_bound = edge;
                has_left = true;
            }
            else
            {
                // Choose the more restrictive (later entry)
                // The intersection of line with left_bound should be "before"
                // the intersection with edge, or vice versa.
                // Compare: at the current left_bound crossing, is the new edge
                // still positive (meaning new edge is later)?
                auto cross_pt = ipg::intersect(line, left_bound);
                auto c = exact_classify(cross_pt, edge);
                if (c > 0)
                {
                    // Current left_bound crossing is still outside new edge
                    // → new edge is more restrictive
                    left_bound = edge;
                }
            }
        }
        else
        {
            // Line exits negative half-space: right boundary
            if (!has_right)
            {
                right_bound = edge;
                has_right = true;
            }
            else
            {
                auto cross_pt = ipg::intersect(line, right_bound);
                auto c = exact_classify(cross_pt, edge);
                if (c > 0)
                {
                    right_bound = edge;
                }
            }
        }
    }

    // Compute endpoints
    if (has_left && has_right)
    {
        // Both boundaries: segment from left crossing to right crossing
        result.v0 = ipg::intersect(line, left_bound);
        result.v1 = ipg::intersect(line, right_bound);

        // Check that left crossing is inside right boundary and vice versa
        auto c0 = exact_classify(result.v0, right_bound);
        if (c0 > 0)
            return {.is_valid = false}; // Empty intersection

        result.is_valid = true;
    }
    else if (has_left && !has_right)
    {
        // Semi-infinite: should not happen for bounded polygon
        return result;
    }
    else if (!has_left && has_right)
    {
        return result;
    }
    else
    {
        // No boundaries: line lies entirely inside polygon
        // This shouldn't happen for a bounded polygon with at least 3 edges
        return result;
    }

    return result;
}

// Compute intersection between two convex polygons
inline PairwiseIntersection intersect_polygons(ConvexPolygon const& t,
                                                ConvexPolygon const& t_prime,
                                                int t_prime_local_idx)
{
    PairwiseIntersection result;

    // Check if supporting planes are parallel
    if (ipg::are_parallel<geometry_t>(t.support, t_prime.support))
    {
        // Parallel planes: check if coplanar (same plane)
        // Two parallel planes are coplanar iff any point of one is on the other
        // Use a vertex of t_prime and classify against t.support
        if (t_prime.vertex_count() == 0) return result;

        auto v0 = t_prime.vertex(0);
        auto c = exact_classify(v0, t.support);
        if (c != 0)
        {
            // Different planes, no intersection
            result.type = PairwiseIntersection::Type::None;
            return result;
        }

        // Coplanar: C4 overlap case
        // Check if polygons actually overlap by testing edge containment
        // (For now, mark as overlap and let BSP handle it)
        result.type = PairwiseIntersection::Type::Overlap;
        result.overlap.other_polygon_idx = t_prime_local_idx;
        result.overlap.other_edges = t_prime.edges;
        result.overlap.other_support = t_prime.support;
        return result;
    }

    // Non-parallel: compute intersection line
    auto int_line = ipg::intersect<geometry_t>(t.support, t_prime.support);

    // Clip line against polygon t's edges
    // We track left and right boundaries
    int n_t = t.vertex_count();
    bool has_left_t = false, has_right_t = false;
    plane_t left_t, right_t;

    for (int i = 0; i < n_t; i++)
    {
        auto const& edge = t.edges[i];
        if (ipg::are_parallel<geometry_t>(edge, int_line))
        {
            auto p = ipg::any_point<geometry_t>(int_line);
            if (exact_classify(p, edge) > 0)
            {
                result.type = PairwiseIntersection::Type::None;
                return result;
            }
            continue;
        }

        auto dir = int_line.direction();
        auto en = edge.normal();
        static constexpr int dot_bits = 3 * geometry_t::bits_normal + 3;
        auto dot = ipg::mul<dot_bits>(dir.x, en.x) +
                   ipg::mul<dot_bits>(dir.y, en.y) +
                   ipg::mul<dot_bits>(dir.z, en.z);

        if (tg::sign(dot) > 0)
        {
            if (!has_left_t) { left_t = edge; has_left_t = true; }
            else
            {
                auto cp = ipg::intersect(int_line, left_t);
                if (exact_classify(cp, edge) > 0)
                    left_t = edge;
            }
        }
        else
        {
            if (!has_right_t) { right_t = edge; has_right_t = true; }
            else
            {
                auto cp = ipg::intersect(int_line, right_t);
                if (exact_classify(cp, edge) > 0)
                    right_t = edge;
            }
        }
    }

    // Clip line against polygon t_prime's edges
    int n_tp = t_prime.vertex_count();
    bool has_left_tp = false, has_right_tp = false;
    plane_t left_tp, right_tp;

    for (int i = 0; i < n_tp; i++)
    {
        auto const& edge = t_prime.edges[i];
        if (ipg::are_parallel<geometry_t>(edge, int_line))
        {
            auto p = ipg::any_point<geometry_t>(int_line);
            if (exact_classify(p, edge) > 0)
            {
                result.type = PairwiseIntersection::Type::None;
                return result;
            }
            continue;
        }

        auto dir = int_line.direction();
        auto en = edge.normal();
        static constexpr int dot_bits = 3 * geometry_t::bits_normal + 3;
        auto dot = ipg::mul<dot_bits>(dir.x, en.x) +
                   ipg::mul<dot_bits>(dir.y, en.y) +
                   ipg::mul<dot_bits>(dir.z, en.z);

        if (tg::sign(dot) > 0)
        {
            if (!has_left_tp) { left_tp = edge; has_left_tp = true; }
            else
            {
                auto cp = ipg::intersect(int_line, left_tp);
                if (exact_classify(cp, edge) > 0)
                    left_tp = edge;
            }
        }
        else
        {
            if (!has_right_tp) { right_tp = edge; has_right_tp = true; }
            else
            {
                auto cp = ipg::intersect(int_line, right_tp);
                if (exact_classify(cp, edge) > 0)
                    right_tp = edge;
            }
        }
    }

    // Merge intervals from both polygons
    // The combined interval is the intersection of both polygon's intervals
    // Each polygon contributes a left and/or right boundary

    // Collect all left boundaries and take the most restrictive
    std::vector<plane_t> lefts, rights;
    if (has_left_t) lefts.push_back(left_t);
    if (has_left_tp) lefts.push_back(left_tp);
    if (has_right_t) rights.push_back(right_t);
    if (has_right_tp) rights.push_back(right_tp);

    if (lefts.empty() || rights.empty())
    {
        // One polygon doesn't bound the line → check if line is fully inside
        // For bounded polygons this shouldn't happen if both have >= 3 edges
        result.type = PairwiseIntersection::Type::None;
        return result;
    }

    // Find most restrictive left boundary
    plane_t best_left = lefts[0];
    for (size_t i = 1; i < lefts.size(); i++)
    {
        auto cp = ipg::intersect(int_line, best_left);
        if (exact_classify(cp, lefts[i]) > 0)
            best_left = lefts[i];
    }

    // Find most restrictive right boundary
    plane_t best_right = rights[0];
    for (size_t i = 1; i < rights.size(); i++)
    {
        auto cp = ipg::intersect(int_line, best_right);
        if (exact_classify(cp, rights[i]) > 0)
            best_right = rights[i];
    }

    // Compute endpoints
    auto v0 = ipg::intersect(int_line, best_left);
    auto v1 = ipg::intersect(int_line, best_right);

    // Check that the interval is non-empty: v0 must be inside right boundary
    auto c_check = exact_classify(v0, best_right);
    if (c_check > 0)
    {
        // Empty intersection
        result.type = PairwiseIntersection::Type::None;
        return result;
    }

    // Check for degenerate case: v0 == v1 (single point)
    if (c_check == 0)
    {
        result.type = PairwiseIntersection::Type::Point;
        return result;
    }

    // Non-degenerate segment
    result.type = PairwiseIntersection::Type::Segment;
    result.segment.v0 = v0;
    result.segment.v1 = v1;
    result.segment.split_plane = t_prime.support;
    result.segment.other_polygon_idx = t_prime_local_idx;
    result.segment.is_valid = true;

    return result;
}

} // namespace ember
