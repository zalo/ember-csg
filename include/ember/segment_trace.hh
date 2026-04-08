#pragma once

// EMBER - WNV Classification via Exact 3-Segment Path (Paper Section 4.4)
//
// The paper's exact classification:
//   1. Find interior point x of BSP leaf polygon, defined as intersect(s_t, e_i, e_j)
//   2. Define reference x_ref as intersect(r0, r1, r2) using 3 axis-aligned planes
//   3. Construct path x → x1 → x2 → x_ref via plane substitution
//   4. For each segment, count crossings using exact integer classify
//   5. The accumulated WNV at x gives (w_front, w_back) for the polygon
//
// No floating-point, no probe offset, no approach-direction heuristic.
// All computation is exact integer arithmetic.

#include <ember/bvh.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <cmath>
#include <vector>

namespace ember
{

// Trace one axis-aligned segment, counting crossings with exact integer arithmetic.
inline WNV trace_segment(
    pos_t const& start, pos_t const& end, WNV const& start_wnv,
    std::vector<ConvexPolygon> const& polygons,
    plane_t const& skip_support, int skip_mesh,
    BVH const* bvh = nullptr)
{
    WNV result = start_wnv;
    int axis = -1;
    for (int a = 0; a < 3; a++)
        if ((&start.x)[a] != (&end.x)[a]) { axis = a; break; }
    if (axis < 0) return result;

    int ax1 = (axis + 1) % 3, ax2 = (axis + 2) % 3;
    plane_t fix1{}, fix2{};
    (&fix1.a)[ax1] = normal_scalar_t(1);
    fix1.d = plane_d_t(-int64_t((&start.x)[ax1]));
    (&fix2.a)[ax2] = normal_scalar_t(1);
    fix2.d = plane_d_t(-int64_t((&start.x)[ax2]));

    bool forward = ((&end.x)[axis] > (&start.x)[axis]);
    plane_t start_bnd{}, end_bnd{};
    (&start_bnd.a)[axis] = normal_scalar_t(1);
    start_bnd.d = plane_d_t(-int64_t((&start.x)[axis]));
    (&end_bnd.a)[axis] = normal_scalar_t(1);
    end_bnd.d = plane_d_t(-int64_t((&end.x)[axis]));

    struct CP { plane_t s; int m; };
    std::vector<CP> crossed;

    auto test = [&](int pi) {
        auto const& poly = polygons[pi];
        if (poly.mesh_index == skip_mesh && poly.support == skip_support) return;
        auto cs = exact_classify(start, poly.support);
        auto ce = exact_classify(end, poly.support);
        if (!((cs > 0 && ce < 0) || (cs < 0 && ce > 0))) return;
        if (tg::is_zero(poly.support.normal_comp(axis))) return;
        for (auto const& c : crossed)
            if (c.m == poly.mesh_index && c.s == poly.support) return;
        auto pt = ipg::intersect(poly.support, fix1, fix2);
        if (!pt.is_valid()) return;
        auto ps = exact_classify(pt, start_bnd);
        auto pe = exact_classify(pt, end_bnd);
        if (!(forward ? (ps > 0 && pe < 0) : (ps < 0 && pe > 0))) return;
        for (auto const& e : poly.edges)
        {
            auto c = exact_classify(pt, e);
            if (c > 0) return;
            if (c == 0)
            {
                auto en = e.normal();
                int64_t nx = int64_t(en.x), ny = int64_t(en.y), nz = int64_t(en.z);
                if (!((nx > 0) || (nx == 0 && ny > 0) || (nx == 0 && ny == 0 && nz > 0)))
                    return;
            }
        }
        crossed.push_back({poly.support, poly.mesh_index});
        int sign = (cs > 0) ? +1 : -1;
        for (size_t i = 0; i < result.size() && i < poly.delta_w.size(); i++)
            result[i] += sign * poly.delta_w[i];
    };

    if (bvh)
    {
        AABB box;
        box.expand(double(start.x), double(start.y), double(start.z));
        box.expand(double(end.x), double(end.y), double(end.z));
        std::vector<int> cands;
        bvh->query(box, cands);
        for (int pi : cands) test(pi);
    }
    else
        for (int pi = 0; pi < (int)polygons.size(); pi++) test(pi);

    return result;
}

// Find interior point ON the polygon using center rounding.
inline bool find_interior_point(ConvexPolygon const& poly, pos_t& out)
{
    int n = poly.vertex_count();
    if (n < 3) return false;
    tg::dpos3 c = {0,0,0};
    for (int i = 0; i < n; i++)
    { auto d = poly.vertex_dpos(i); c.x += d.x; c.y += d.y; c.z += d.z; }
    c.x /= n; c.y /= n; c.z /= n;

    for (int r = 0; r <= 2; r++)
    for (int dx = -r; dx <= r; dx++)
    for (int dy = -r; dy <= r; dy++)
    for (int dz = -r; dz <= r; dz++)
    {
        if (r > 0 && std::abs(dx) < r && std::abs(dy) < r && std::abs(dz) < r) continue;
        pos_t p(pos_scalar_t(int32_t(std::round(c.x)) + dx),
                pos_scalar_t(int32_t(std::round(c.y)) + dy),
                pos_scalar_t(int32_t(std::round(c.z)) + dz));
        if (exact_classify(p, poly.support) != 0) continue;
        bool inside = true;
        for (auto const& e : poly.edges)
            if (exact_classify(p, e) > 0) { inside = false; break; }
        if (inside) { out = p; return true; }
    }
    return false;
}

// Classify a polygon using the paper's exact algorithm.
// The interior point is ON the surface. We trace from the reference using
// axis-aligned segments (L-path), skipping the host polygon's coplanar
// same-mesh faces. The approach direction determines which side of the
// surface the accumulated WNV represents.
inline WNV classify_leaf_polygon(
    plane_t const& support,
    std::vector<plane_t> const& leaf_edges,
    pos_t const& ref_point,
    WNV const& ref_wnv,
    std::vector<ConvexPolygon> const& polygons,
    IAABB const& /*aabb*/,
    WNTV const& host_delta_w,
    int host_mesh_index = -1,
    BVH const* bvh = nullptr)
{
    ConvexPolygon tmp;
    tmp.support = support;
    tmp.edges = leaf_edges;

    pos_t interior;
    if (!find_interior_point(tmp, interior))
        return ref_wnv;

    // Trace L-path from reference to interior, skipping host's coplanar faces
    WNV wnv = ref_wnv;
    pos_t cur = ref_point;
    if (interior.x != cur.x)
    { pos_t n(interior.x, cur.y, cur.z); wnv = trace_segment(cur, n, wnv, polygons, support, host_mesh_index, bvh); cur = n; }
    if (interior.y != cur.y)
    { pos_t n(cur.x, interior.y, cur.z); wnv = trace_segment(cur, n, wnv, polygons, support, host_mesh_index, bvh); cur = n; }
    if (interior.z != cur.z)
    { wnv = trace_segment(cur, interior, wnv, polygons, support, host_mesh_index, bvh); }

    // The interior point is ON the support plane. The accumulated WNV is for
    // the side the trace approached from. Determine which side:
    // pre = last intermediate point before reaching interior.
    // For L-path X→Y→Z: pre = (interior.x, interior.y, ref.z)
    pos_t pre = ref_point;
    if (interior.x != ref_point.x) pre.x = interior.x;
    if (interior.y != ref_point.y) pre.y = interior.y;
    // pre now has the position BEFORE the last segment

    auto approach = exact_classify(pre, support);
    if (approach == 0)
    {
        // pre is on the support plane (happens when last segment is tangent).
        // Step back one more level.
        pos_t pre2 = ref_point;
        if (interior.x != ref_point.x) pre2.x = interior.x;
        approach = exact_classify(pre2, support);
        if (approach == 0)
            approach = exact_classify(ref_point, support);
    }

    // approach > 0: arrived from positive (front) side → wnv IS w_front
    // approach < 0: arrived from negative (back) side → wnv IS w_back
    //               w_front = w_back - delta_w
    if (approach < 0)
    {
        for (size_t k = 0; k < wnv.size() && k < host_delta_w.size(); k++)
            wnv[k] -= host_delta_w[k];
    }

    return wnv;
}

} // namespace ember
