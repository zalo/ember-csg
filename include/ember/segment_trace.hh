#pragma once

// EMBER - WNV Classification via robust multi-ray casting + exact segment trace.
//
// Two modes:
//   1. point_in_meshes_robust: 7-ray majority vote, BVH-accelerated (for leaf classification)
//   2. trace_axis_segment: exact integer segment trace (for reference propagation)

#include <ember/bvh.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/exact_classify.hh>
#include <ember/winding.hh>

#include <cmath>
#include <vector>

namespace ember
{

// ---------------------------------------------------------------------------
// Axis-aligned plane + exact segment trace (for reference propagation)
// ---------------------------------------------------------------------------

inline plane_t axis_plane(int axis, pos_scalar_t val)
{
    plane_t p;
    p.a = normal_scalar_t(axis == 0 ? 1 : 0);
    p.b = normal_scalar_t(axis == 1 ? 1 : 0);
    p.c = normal_scalar_t(axis == 2 ? 1 : 0);
    p.d = plane_d_t(-int64_t(val));
    return p;
}

// Exact axis-aligned segment trace. Accumulates WNV transitions for all
// polygon crossings along the segment.
inline WNV trace_axis_segment(
    pos_t const& start, pos_t const& end, int axis,
    WNV const& start_wnv,
    std::vector<ConvexPolygon> const& polygons,
    bool& valid)
{
    valid = true;
    WNV wnv = start_wnv;

    int ax1 = (axis + 1) % 3;
    int ax2 = (axis + 2) % 3;
    plane_t lp0 = axis_plane(ax1, (&start.x)[ax1]);
    plane_t lp1 = axis_plane(ax2, (&start.x)[ax2]);

    int dir_sign = (int32_t((&start.x)[axis]) < int32_t((&end.x)[axis])) ? 1 : -1;
    plane_t bound_start = axis_plane(axis, (&start.x)[axis]);
    plane_t bound_end   = axis_plane(axis, (&end.x)[axis]);

    for (auto const& poly : polygons)
    {
        auto sn = poly.support.normal();
        if (int64_t((&sn.x)[axis]) == 0) continue;

        auto pt = ipg::intersect(lp0, lp1, poly.support);
        if (!pt.is_valid()) continue;

        auto cs = exact_classify(pt, bound_start);
        auto ce = exact_classify(pt, bound_end);

        if (dir_sign > 0) { if (cs <= 0 || ce >= 0) continue; }
        else              { if (cs >= 0 || ce <= 0) continue; }

        bool inside = true;
        bool on_edge = false;
        for (auto const& edge : poly.edges)
        {
            auto c = exact_classify(pt, edge);
            if (c > 0) { inside = false; break; }
            if (c == 0) on_edge = true;
        }
        if (!inside) continue;
        if (on_edge) { valid = false; return wnv; }

        int cross_sign = (int64_t((&sn.x)[axis]) > 0 ? 1 : -1) * (-dir_sign);
        for (size_t k = 0; k < wnv.size() && k < poly.delta_w.size(); k++)
            wnv[k] += cross_sign * poly.delta_w[k];
    }
    return wnv;
}

// trace_segment: exact L-path trace for reference propagation.
inline WNV trace_segment(pos_t const& start, pos_t const& end, WNV const& wnv,
    std::vector<ConvexPolygon> const& polys, plane_t const&, int,
    BVH const* = nullptr)
{
    static const int orderings[][3] = {{0,1,2}, {1,2,0}, {2,0,1}};
    for (auto& ord : orderings)
    {
        bool valid = true;
        WNV attempt = wnv;
        pos_t cur = start;
        for (int i = 0; i < 3 && valid; i++)
        {
            int ax = ord[i];
            if ((&cur.x)[ax] != (&end.x)[ax])
            {
                pos_t next = cur;
                (&next.x)[ax] = (&end.x)[ax];
                attempt = trace_axis_segment(cur, next, ax, attempt, polys, valid);
                cur = next;
            }
        }
        if (valid) return attempt;
    }
    return wnv; // fallback
}

// ---------------------------------------------------------------------------
// 7-ray robust point-in-mesh (for leaf polygon classification)
// ---------------------------------------------------------------------------

inline double ray_tri_hit(double ox,double oy,double oz,
                          double dx,double dy,double dz,
                          double v0x,double v0y,double v0z,
                          double v1x,double v1y,double v1z,
                          double v2x,double v2y,double v2z)
{
    double e1x=v1x-v0x,e1y=v1y-v0y,e1z=v1z-v0z;
    double e2x=v2x-v0x,e2y=v2y-v0y,e2z=v2z-v0z;
    double hx=dy*e2z-dz*e2y,hy=dz*e2x-dx*e2z,hz=dx*e2y-dy*e2x;
    double a=e1x*hx+e1y*hy+e1z*hz;
    if(std::abs(a)<1e-30) return -1;
    double f=1.0/a;
    double sx=ox-v0x,sy=oy-v0y,sz=oz-v0z;
    double u=f*(sx*hx+sy*hy+sz*hz);
    if(u<-1e-8||u>1+1e-8) return -1;
    double qx=sy*e1z-sz*e1y,qy=sz*e1x-sx*e1z,qz=sx*e1y-sy*e1x;
    double v=f*(dx*qx+dy*qy+dz*qz);
    if(v<-1e-8||u+v>1+1e-8) return -1;
    return f*(e2x*qx+e2y*qy+e2z*qz);
}

// 7-ray majority vote, BVH-accelerated
inline WNV point_in_meshes_robust(
    double px, double py, double pz,
    std::vector<ConvexPolygon> const& polygons, int num_meshes,
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
        auto test = [&](int pi) {
            auto const& poly = polygons[pi];
            if (poly.mesh_index < 0 || poly.mesh_index >= num_meshes) return;
            int nv = poly.vertex_count();
            if (nv < 3) return;
            auto v0 = poly.vertex_dpos(0);
            for (int k = 1; k < nv-1; k++) {
                auto v1 = poly.vertex_dpos(k);
                auto v2 = poly.vertex_dpos(k+1);
                double t = ray_tri_hit(px,py,pz, d[0],d[1],d[2],
                    v0.x,v0.y,v0.z, v1.x,v1.y,v1.z, v2.x,v2.y,v2.z);
                if (t > 1e-6) crossings[poly.mesh_index]++;
            }
        };
        if (bvh)
            bvh->raycast(px,py,pz, d[0],d[1],d[2], [&](int pi){ test(pi); });
        else
            for (int pi = 0; pi < (int)polygons.size(); pi++) test(pi);
        for (int m = 0; m < num_meshes; m++)
            if (crossings[m] & 1) votes[m]++;
    }
    WNV result(num_meshes, 0);
    for (int m = 0; m < num_meshes; m++)
        result[m] = (votes[m] > 3) ? 1 : 0;
    return result;
}

// Probe: center + 10% edge offset along dominant normal
inline bool find_probe_point(ConvexPolygon const& poly, pos_t& out_probe, int& out_side)
{
    int n = poly.vertex_count();
    if (n < 3) return false;
    tg::dpos3 center = {0,0,0};
    for (int i = 0; i < n; i++)
    { auto dp = poly.vertex_dpos(i); center.x+=dp.x; center.y+=dp.y; center.z+=dp.z; }
    center.x/=n; center.y/=n; center.z/=n;

    auto sn = poly.support.normal();
    auto aa = ipg::abs(int64_t(sn.x)), ab = ipg::abs(int64_t(sn.y)), ac = ipg::abs(int64_t(sn.z));
    int dom = 0;
    if (ab > aa && ab > ac) dom = 1;
    else if (ac > aa && ac > ab) dom = 2;
    int dir = (int64_t((&sn.x)[dom]) > 0) ? 1 : -1;

    double edge_len = 0;
    for (int i = 0; i < n; i++) {
        auto vi = poly.vertex_dpos(i); auto vj = poly.vertex_dpos((i+1)%n);
        double dx=vi.x-vj.x,dy=vi.y-vj.y,dz=vi.z-vj.z;
        edge_len = std::max(edge_len, std::sqrt(dx*dx+dy*dy+dz*dz));
    }
    int offset = std::max(int32_t(edge_len * 0.1), int32_t(2));

    pos_t probe(pos_scalar_t(int32_t(std::round(center.x))),
                pos_scalar_t(int32_t(std::round(center.y))),
                pos_scalar_t(int32_t(std::round(center.z))));
    (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir * offset);
    auto side = exact_classify(probe, poly.support);
    if (side == 0) { (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir*offset); side = exact_classify(probe, poly.support); }
    if (side == 0) return false;
    out_probe = probe; out_side = side;
    return true;
}

inline bool find_interior_point(ConvexPolygon const& poly, pos_t& out)
{
    int n = poly.vertex_count();
    if (n < 3) return false;
    tg::dpos3 c = {0,0,0};
    for (int i = 0; i < n; i++) { auto d = poly.vertex_dpos(i); c.x+=d.x; c.y+=d.y; c.z+=d.z; }
    c.x/=n; c.y/=n; c.z/=n;
    for (int r = 0; r <= 2; r++)
    for (int dx = -r; dx <= r; dx++)
    for (int dy = -r; dy <= r; dy++)
    for (int dz = -r; dz <= r; dz++) {
        if (r>0&&std::abs(dx)<r&&std::abs(dy)<r&&std::abs(dz)<r) continue;
        pos_t p(pos_scalar_t(int32_t(std::round(c.x))+dx),pos_scalar_t(int32_t(std::round(c.y))+dy),pos_scalar_t(int32_t(std::round(c.z))+dz));
        if (exact_classify(p,poly.support)!=0) continue;
        bool inside=true;
        for (auto const& e : poly.edges) if (exact_classify(p,e)>0) { inside=false; break; }
        if (inside) { out=p; return true; }
    }
    return false;
}

// Classify a leaf polygon via exact segment tracing from the reference point.
// Finds a probe point off the polygon surface, traces an L-path from ref
// to probe through the local polygon set, accumulating WNV transitions.
inline WNV classify_leaf_polygon(
    plane_t const& support, std::vector<plane_t> const& leaf_edges,
    pos_t const& ref_point, WNV const& ref_wnv,
    std::vector<ConvexPolygon> const& polygons, IAABB const& aabb,
    WNTV const& host_delta_w, int /*host_mesh*/ = -1, BVH const* = nullptr)
{
    ConvexPolygon tmp;
    tmp.support = support;
    tmp.edges = leaf_edges;

    pos_t probe;
    int probe_side;
    if (!find_probe_point(tmp, probe, probe_side))
        return ref_wnv;

    // Clamp probe to AABB so the L-path stays within the local polygon region
    for (int a = 0; a < 3; a++)
    {
        auto& coord = (&probe.x)[a];
        auto bmin = (&aabb.min.x)[a] + pos_scalar_t(1);
        auto bmax = (&aabb.max.x)[a] - pos_scalar_t(1);
        if (bmin < bmax)
        {
            if (coord < bmin) coord = bmin;
            else if (coord > bmax) coord = bmax;
        }
    }
    // Recompute side after potential clamping
    probe_side = exact_classify(probe, support);
    if (probe_side == 0)
        return ref_wnv;

    // Trace L-path from reference to probe, trying multiple axis orderings
    static const int orderings[][3] = {{0,1,2}, {1,2,0}, {2,0,1}};
    for (auto& ord : orderings)
    {
        bool valid = true;
        WNV wnv = ref_wnv;
        pos_t cur = ref_point;
        for (int i = 0; i < 3 && valid; i++)
        {
            int ax = ord[i];
            if ((&cur.x)[ax] != (&probe.x)[ax])
            {
                pos_t next = cur;
                (&next.x)[ax] = (&probe.x)[ax];
                wnv = trace_axis_segment(cur, next, ax, wnv, polygons, valid);
                cur = next;
            }
        }
        if (valid)
        {
            // Side correction
            if (probe_side < 0)
                for (size_t k = 0; k < wnv.size() && k < host_delta_w.size(); k++)
                    wnv[k] -= host_delta_w[k];
            return wnv;
        }
    }

    return ref_wnv; // All paths hit edges — fallback
}

} // namespace ember
