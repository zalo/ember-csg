#pragma once

// EMBER - WNV Classification via robust multi-ray casting.
//
// Matches the aa1b8b6 "manifold-producing" state, with BVH acceleration.
//
// Algorithm:
//   1. Find probe point: polygon center + 10% edge-length offset along normal
//   2. Cast 7 rays in random directions, count parity per mesh (majority vote)
//   3. NO coplanar skip (probe is off-surface, all crossings are real)
//   4. Side correction: if probe on back (negative) side, subtract delta_w

#include <ember/bvh.hh>
#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <cmath>
#include <vector>

namespace ember
{

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

// 7-ray majority vote, NO coplanar skip, BVH-accelerated
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

// Classify using aa1b8b6 logic: probe + robust multi-ray + side correction
inline WNV classify_leaf_polygon(
    plane_t const& support, std::vector<plane_t> const& leaf_edges,
    pos_t const& /*ref*/, WNV const& ref_wnv,
    std::vector<ConvexPolygon> const& polygons, IAABB const& /*aabb*/,
    WNTV const& host_delta_w, int /*host_mesh*/ = -1, BVH const* bvh = nullptr)
{
    ConvexPolygon tmp; tmp.support = support; tmp.edges = leaf_edges;
    pos_t probe; int probe_side;
    if (!find_probe_point(tmp, probe, probe_side)) return ref_wnv;
    int nm = (int)ref_wnv.size();
    WNV wnv = point_in_meshes_robust(double(probe.x), double(probe.y), double(probe.z), polygons, nm, bvh);
    if (probe_side > 0) return wnv;
    for (size_t k = 0; k < wnv.size() && k < host_delta_w.size(); k++) wnv[k] -= host_delta_w[k];
    return wnv;
}

// trace_segment for reference propagation
inline WNV trace_segment(pos_t const& start, pos_t const& end, WNV const& wnv,
    std::vector<ConvexPolygon> const& polys, plane_t const&, int, BVH const* bvh = nullptr)
{
    return point_in_meshes_robust(double(end.x),double(end.y),double(end.z), polys, (int)wnv.size(), bvh);
}

} // namespace ember
