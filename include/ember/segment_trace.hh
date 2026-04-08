#pragma once

// EMBER - WNV Classification
//
// Two implementations:
//   1. Robust: BVH-accelerated multi-ray casting (Moller-Trumbore, 7 rays, majority vote)
//   2. Exact:  Paper's 3-segment path via plane substitution (Section 4.4)
//
// Both skip the host polygon's coplanar same-mesh faces.
// The robust method is the reference for validation.

#include <ember/bvh.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <cmath>
#include <vector>

namespace ember
{

// ============================================================================
// ROBUST: Multi-ray casting (BVH-accelerated Moller-Trumbore)
// ============================================================================

inline double ray_tri_hit(double ox, double oy, double oz,
                          double dx, double dy, double dz,
                          tg::dpos3 const& v0, tg::dpos3 const& v1, tg::dpos3 const& v2)
{
    double e1x=v1.x-v0.x, e1y=v1.y-v0.y, e1z=v1.z-v0.z;
    double e2x=v2.x-v0.x, e2y=v2.y-v0.y, e2z=v2.z-v0.z;
    double hx=dy*e2z-dz*e2y, hy=dz*e2x-dx*e2z, hz=dx*e2y-dy*e2x;
    double a=e1x*hx+e1y*hy+e1z*hz;
    if (std::abs(a)<1e-30) return -1;
    double f=1.0/a;
    double sx=ox-v0.x, sy=oy-v0.y, sz=oz-v0.z;
    double u=f*(sx*hx+sy*hy+sz*hz);
    if (u<-1e-8||u>1+1e-8) return -1;
    double qx=sy*e1z-sz*e1y, qy=sz*e1x-sx*e1z, qz=sx*e1y-sy*e1x;
    double v=f*(dx*qx+dy*qy+dz*qz);
    if (v<-1e-8||u+v>1+1e-8) return -1;
    return f*(e2x*qx+e2y*qy+e2z*qz);
}

inline WNV robust_classify(
    double px, double py, double pz,
    std::vector<ConvexPolygon> const& polygons, int num_meshes,
    plane_t const& skip_support, int skip_mesh,
    BVH const* bvh = nullptr)
{
    static const double dirs[][3] = {
        { 0.85065, 0.52573, 0.03532}, {-0.38912, 0.92131, 0.01234},
        { 0.12340,-0.45670, 0.88092}, {-0.70711,-0.70711, 0.01002},
        { 0.57735, 0.57735, 0.57735}, {-0.23450, 0.12340,-0.96421},
        { 0.95106, 0.03123,-0.30734},
    };

    std::vector<int> votes(num_meshes, 0);
    for (auto& d : dirs)
    {
        std::vector<int> crossings(num_meshes, 0);
        auto test_poly = [&](int pi) {
            auto const& poly = polygons[pi];
            if (poly.mesh_index == skip_mesh && poly.support == skip_support) return;
            if (poly.mesh_index < 0 || poly.mesh_index >= num_meshes) return;
            int nv = poly.vertex_count();
            if (nv < 3) return;
            auto v0 = poly.vertex_dpos(0);
            for (int k = 1; k < nv-1; k++)
            {
                auto v1 = poly.vertex_dpos(k);
                auto v2 = poly.vertex_dpos(k+1);
                double t = ray_tri_hit(px, py, pz, d[0], d[1], d[2], v0, v1, v2);
                if (t > 1e-6) crossings[poly.mesh_index]++;
            }
        };
        if (bvh)
            bvh->raycast(px, py, pz, d[0], d[1], d[2], [&](int pi) { test_poly(pi); });
        else
            for (int pi = 0; pi < (int)polygons.size(); pi++) test_poly(pi);
        for (int m = 0; m < num_meshes; m++)
            if (crossings[m] & 1) votes[m]++;
    }
    WNV result(num_meshes, 0);
    for (int m = 0; m < num_meshes; m++)
        result[m] = (votes[m] > 3) ? 1 : 0;
    return result;
}

// ============================================================================
// EXACT: 3-segment path via plane substitution (Paper Section 4.4)
// ============================================================================

// Trace a segment on LINE between two point4s bounded by two planes.
// C is strictly between A and B iff:
//   sign(classify(C, bound_a)) == sign(classify(B, bound_a)) AND != 0
//   sign(classify(C, bound_b)) == sign(classify(A, bound_b)) AND != 0
inline WNV trace_line_segment(
    line_t const& line,
    point4_t const& pt_a, point4_t const& pt_b,
    plane_t const& bound_a, plane_t const& bound_b,
    WNV const& start_wnv,
    std::vector<ConvexPolygon> const& polygons,
    plane_t const& skip_support, int skip_mesh,
    BVH const* bvh = nullptr)
{
    WNV result = start_wnv;
    auto b_vs_ba = exact_classify(pt_b, bound_a);
    auto a_vs_bb = exact_classify(pt_a, bound_b);
    if (b_vs_ba == 0 || a_vs_bb == 0) return result;

    struct CP { plane_t s; int m; };
    std::vector<CP> crossed;

    auto test = [&](int pi) {
        auto const& poly = polygons[pi];
        if (poly.mesh_index == skip_mesh && poly.support == skip_support) return;
        if (ipg::are_parallel<geometry_t>(poly.support, line)) return;
        for (auto const& c : crossed)
            if (c.m == poly.mesh_index && c.s == poly.support) return;

        auto cross_pt = ipg::intersect(line, poly.support);
        if (!cross_pt.is_valid()) return;

        auto c_ba = exact_classify(cross_pt, bound_a);
        auto c_bb = exact_classify(cross_pt, bound_b);
        if (c_ba == 0 || c_bb == 0) return;
        if ((c_ba > 0) != (b_vs_ba > 0)) return;
        if ((c_bb > 0) != (a_vs_bb > 0)) return;

        for (auto const& e : poly.edges)
        {
            auto c = exact_classify(cross_pt, e);
            if (c > 0) return;
            if (c == 0)
            {
                auto en = e.normal();
                int64_t nx=int64_t(en.x), ny=int64_t(en.y), nz=int64_t(en.z);
                if (!((nx>0)||(nx==0&&ny>0)||(nx==0&&ny==0&&nz>0))) return;
            }
        }

        auto sa = exact_classify(pt_a, poly.support);
        int sign = (sa > 0) ? +1 : (sa < 0) ? -1 : 0;
        if (sign == 0) return;

        crossed.push_back({poly.support, poly.mesh_index});
        for (size_t i = 0; i < result.size() && i < poly.delta_w.size(); i++)
            result[i] += sign * poly.delta_w[i];
    };

    if (bvh)
    {
        auto da = ipg::to_dpos3(pt_a);
        auto db = ipg::to_dpos3(pt_b);
        AABB box;
        box.expand(da.x, da.y, da.z);
        box.expand(db.x, db.y, db.z);
        for (int a = 0; a < 3; a++) { box.min[a] -= 100; box.max[a] += 100; }
        std::vector<int> cands;
        bvh->query(box, cands);
        for (int pi : cands) test(pi);
    }
    else
        for (int pi = 0; pi < (int)polygons.size(); pi++) test(pi);

    return result;
}

// Axis-aligned segment trace (for reference propagation)
inline WNV trace_segment(
    pos_t const& start, pos_t const& end, WNV const& start_wnv,
    std::vector<ConvexPolygon> const& polygons,
    plane_t const& skip_support, int skip_mesh,
    BVH const* bvh = nullptr)
{
    int axis = -1;
    for (int a = 0; a < 3; a++)
        if ((&start.x)[a] != (&end.x)[a]) { axis = a; break; }
    if (axis < 0) return start_wnv;

    int ax1=(axis+1)%3, ax2=(axis+2)%3;
    plane_t fix1{}, fix2{};
    (&fix1.a)[ax1] = normal_scalar_t(1);
    fix1.d = plane_d_t(-int64_t((&start.x)[ax1]));
    (&fix2.a)[ax2] = normal_scalar_t(1);
    fix2.d = plane_d_t(-int64_t((&start.x)[ax2]));

    plane_t sb{}, eb{};
    (&sb.a)[axis] = normal_scalar_t(1);
    sb.d = plane_d_t(-int64_t((&start.x)[axis]));
    (&eb.a)[axis] = normal_scalar_t(1);
    eb.d = plane_d_t(-int64_t((&end.x)[axis]));

    auto line = ipg::intersect<geometry_t>(fix1, fix2);
    point4_t hs(tg::ipos3(int(start.x), int(start.y), int(start.z)));
    point4_t he(tg::ipos3(int(end.x), int(end.y), int(end.z)));

    return trace_line_segment(line, hs, he, sb, eb, start_wnv, polygons, skip_support, skip_mesh, bvh);
}

// ============================================================================
// INTERIOR POINT FINDING
// ============================================================================

inline bool find_interior_point(ConvexPolygon const& poly, pos_t& out,
                                 int& out_edge_i, int& out_edge_j)
{
    int n = poly.vertex_count();
    if (n < 3) return false;
    tg::dpos3 c = {0,0,0};
    for (int i = 0; i < n; i++)
    { auto d = poly.vertex_dpos(i); c.x+=d.x; c.y+=d.y; c.z+=d.z; }
    c.x/=n; c.y/=n; c.z/=n;

    for (int r = 0; r <= 2; r++)
    for (int dx = -r; dx <= r; dx++)
    for (int dy = -r; dy <= r; dy++)
    for (int dz = -r; dz <= r; dz++)
    {
        if (r>0 && std::abs(dx)<r && std::abs(dy)<r && std::abs(dz)<r) continue;
        pos_t p(pos_scalar_t(int32_t(std::round(c.x))+dx),
                pos_scalar_t(int32_t(std::round(c.y))+dy),
                pos_scalar_t(int32_t(std::round(c.z))+dz));
        if (exact_classify(p, poly.support) != 0) continue;
        bool inside = true;
        for (auto const& e : poly.edges)
            if (exact_classify(p, e) > 0) { inside = false; break; }
        if (inside)
        {
            out = p;
            double best = 1e30;
            out_edge_i = 0; out_edge_j = 1;
            for (int k = 0; k < n; k++)
            {
                auto vk = poly.vertex_dpos(k);
                double d2 = (vk.x-double(p.x))*(vk.x-double(p.x)) +
                            (vk.y-double(p.y))*(vk.y-double(p.y)) +
                            (vk.z-double(p.z))*(vk.z-double(p.z));
                if (d2 < best) { best=d2; out_edge_i=k; out_edge_j=(k+1)%n; }
            }
            return true;
        }
    }
    return false;
}

inline bool find_interior_point(ConvexPolygon const& poly, pos_t& out)
{ int ei, ej; return find_interior_point(poly, out, ei, ej); }

// ============================================================================
// MAIN CLASSIFY: tries exact 3-segment, falls back to robust multi-ray
// ============================================================================

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
    int edge_i, edge_j;
    if (!find_interior_point(tmp, interior, edge_i, edge_j))
        return ref_wnv;

    int num_meshes = (int)ref_wnv.size();

    // --- PRIMARY: robust multi-ray at an off-surface probe point ---
    // Offset from surface along the dominant normal axis by ~10% of edge length.
    // This keeps the probe close to the polygon (correct region) while being
    // far enough from the surface to avoid edge-clipping ambiguities.
    auto sn = support.normal();
    auto abs_a = ipg::abs(int64_t(sn.x)), abs_b = ipg::abs(int64_t(sn.y)), abs_c = ipg::abs(int64_t(sn.z));
    int dom = 0;
    if (abs_b > abs_a && abs_b > abs_c) dom = 1;
    else if (abs_c > abs_a && abs_c > abs_b) dom = 2;
    int dir = (int64_t((&sn.x)[dom]) > 0) ? 1 : -1;

    double max_edge = 0;
    int ne = tmp.vertex_count();
    for (int i = 0; i < ne; i++)
    {
        auto vi = tmp.vertex_dpos(i);
        auto vj = tmp.vertex_dpos((i+1) % ne);
        double dx = vi.x-vj.x, dy = vi.y-vj.y, dz = vi.z-vj.z;
        max_edge = std::max(max_edge, std::sqrt(dx*dx+dy*dy+dz*dz));
    }
    int offset = std::max(int32_t(max_edge * 0.1), int32_t(2));

    pos_t probe = interior;
    (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir * offset);
    auto side = exact_classify(probe, support);
    if (side == 0) { (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir * offset); side = exact_classify(probe, support); }
    if (side == 0) return ref_wnv;

    double px = double(probe.x), py = double(probe.y), pz = double(probe.z);
    WNV wnv = robust_classify(px, py, pz, polygons, num_meshes,
                                support, host_mesh_index, bvh);

    if (side < 0)
        for (size_t k = 0; k < wnv.size() && k < host_delta_w.size(); k++)
            wnv[k] -= host_delta_w[k];

    return wnv;
}

} // namespace ember
