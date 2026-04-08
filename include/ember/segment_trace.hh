#pragma once

// EMBER - WNV Classification via robust double-precision ray casting.
//
// Uses Moller-Trumbore ray-triangle intersection with multiple ray directions
// and majority voting for robust inside/outside determination.
// The exact integer arithmetic is used for everything else (plane intersections,
// polygon clipping, BSP construction). Only the final classification step
// uses double precision, matching the paper's note that classification
// can use fast approximate methods when the geometry is exact.

#include <ember/polygon.hh>
#include <ember/types.hh>
#include <ember/winding.hh>

#include <cmath>
#include <cstdio>
#include <vector>

namespace ember
{

// Moller-Trumbore ray-triangle intersection (double precision)
inline double ray_tri_hit(double ox, double oy, double oz,
                          double dx, double dy, double dz,
                          double v0x, double v0y, double v0z,
                          double v1x, double v1y, double v1z,
                          double v2x, double v2y, double v2z)
{
    double e1x = v1x-v0x, e1y = v1y-v0y, e1z = v1z-v0z;
    double e2x = v2x-v0x, e2y = v2y-v0y, e2z = v2z-v0z;
    double hx = dy*e2z - dz*e2y, hy = dz*e2x - dx*e2z, hz = dx*e2y - dy*e2x;
    double a = e1x*hx + e1y*hy + e1z*hz;
    if (std::abs(a) < 1e-30) return -1;
    double f = 1.0/a;
    double sx = ox-v0x, sy = oy-v0y, sz = oz-v0z;
    double u = f * (sx*hx + sy*hy + sz*hz);
    if (u < -1e-8 || u > 1.0+1e-8) return -1;
    double qx = sy*e1z - sz*e1y, qy = sz*e1x - sx*e1z, qz = sx*e1y - sy*e1x;
    double v = f * (dx*qx + dy*qy + dz*qz);
    if (v < -1e-8 || u+v > 1.0+1e-8) return -1;
    return f * (e2x*qx + e2y*qy + e2z*qz);
}

// Cast one ray from point, count parity of crossings per mesh.
// Uses the original triangle vertex positions (from vertex_dpos) for accuracy.
inline void cast_ray_parity(
    double px, double py, double pz,
    double dx, double dy, double dz,
    std::vector<ConvexPolygon> const& polygons,
    int num_meshes,
    std::vector<int>& crossings) // crossings[mesh] += hit count
{
    for (auto const& poly : polygons)
    {
        if (poly.mesh_index < 0 || poly.mesh_index >= num_meshes) continue;
        int nv = poly.vertex_count();
        if (nv < 3) continue;

        // Fan-triangulate and test each sub-triangle
        auto v0 = poly.vertex_dpos(0);
        for (int k = 1; k < nv - 1; k++)
        {
            auto v1 = poly.vertex_dpos(k);
            auto v2 = poly.vertex_dpos(k + 1);
            double t = ray_tri_hit(px, py, pz, dx, dy, dz,
                v0.x, v0.y, v0.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
            if (t > 1e-6)
                crossings[poly.mesh_index]++;
        }
    }
}

// Robust point-in-meshes test: cast 7 rays in different directions,
// take majority vote per mesh for inside/outside.
inline WNV point_in_meshes_robust(
    double px, double py, double pz,
    std::vector<ConvexPolygon> const& polygons,
    int num_meshes)
{
    static const double dirs[][3] = {
        { 0.85065, 0.52573, 0.03532},
        {-0.38912, 0.92131, 0.01234},
        { 0.12340,-0.45670, 0.88092},
        {-0.70711,-0.70711, 0.01002},
        { 0.57735, 0.57735, 0.57735},
        {-0.23450, 0.12340,-0.96421},
        { 0.95106, 0.03123,-0.30734},
    };

    std::vector<int> votes(num_meshes, 0);

    for (auto& d : dirs)
    {
        std::vector<int> crossings(num_meshes, 0);
        cast_ray_parity(px, py, pz, d[0], d[1], d[2], polygons, num_meshes, crossings);
        for (int m = 0; m < num_meshes; m++)
            if (crossings[m] & 1) votes[m]++;
    }

    WNV result(num_meshes, 0);
    for (int m = 0; m < num_meshes; m++)
        result[m] = (votes[m] > 3) ? 1 : 0;
    return result;
}

// Find a probe point near a polygon. Uses ConvexPolygon's vertex_dpos
// (binary-search based) for accurate center computation.
inline bool find_probe_point(ConvexPolygon const& poly,
                              pos_t& out_probe,
                              int& out_side)
{
    auto const& support = poly.support;
    int n = poly.vertex_count();
    if (n < 3) return false;

    tg::dpos3 center = {0, 0, 0};
    for (int i = 0; i < n; i++)
    {
        auto dp = poly.vertex_dpos(i);
        center.x += dp.x;
        center.y += dp.y;
        center.z += dp.z;
    }
    center.x /= n;
    center.y /= n;
    center.z /= n;

    auto sn = support.normal();
    auto abs_a = ipg::abs(int64_t(sn.x));
    auto abs_b = ipg::abs(int64_t(sn.y));
    auto abs_c = ipg::abs(int64_t(sn.z));
    int dom = 0;
    if (abs_b > abs_a && abs_b > abs_c) dom = 1;
    else if (abs_c > abs_a && abs_c > abs_b) dom = 2;
    int dir = (int64_t((&sn.x)[dom]) > 0) ? 1 : -1;

    // Estimate triangle size to set a safe probe offset.
    // The probe must be far enough from the surface that rays don't
    // clip through adjacent triangles at shared edges.
    double edge_len = 0;
    for (int i = 0; i < n; i++)
    {
        auto vi = poly.vertex_dpos(i);
        auto vj = poly.vertex_dpos((i+1) % n);
        double dx = vi.x - vj.x, dy = vi.y - vj.y, dz = vi.z - vj.z;
        double len = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (len > edge_len) edge_len = len;
    }
    // Offset by ~10% of the longest edge (enough to clear adjacent triangles)
    int offset = std::max(int32_t(edge_len * 0.1), int32_t(2));

    pos_t probe(pos_scalar_t(int32_t(std::round(center.x))),
                pos_scalar_t(int32_t(std::round(center.y))),
                pos_scalar_t(int32_t(std::round(center.z))));
    (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir * offset);

    auto side = exact_classify(probe, support);
    if (side == 0)
    {
        (&probe.x)[dom] = pos_scalar_t(int32_t((&probe.x)[dom]) + dir * offset);
        side = exact_classify(probe, support);
    }
    if (side == 0) return false;

    out_probe = probe;
    out_side = side;
    return true;
}

// Classify a polygon: find WNV at a probe point near its surface.
inline WNV classify_leaf_polygon(
    plane_t const& support,
    std::vector<plane_t> const& leaf_edges,
    pos_t const& /*ref_point*/,
    WNV const& ref_wnv,
    std::vector<ConvexPolygon> const& polygons,
    IAABB const& aabb,
    WNTV const& host_delta_w)
{
    // Build a temporary ConvexPolygon for probe point finding
    ConvexPolygon tmp_poly;
    tmp_poly.support = support;
    tmp_poly.edges = leaf_edges;

    pos_t probe;
    int probe_side;
    if (!find_probe_point(tmp_poly, probe, probe_side))
        return ref_wnv;

    int num_meshes = static_cast<int>(ref_wnv.size());

    // Convert probe to double for ray casting
    // Use the transform stored in the AABB context... we need original-space coords.
    // For the ray casting, we use integer-space coordinates directly since
    // vertex_dpos also returns integer-space coordinates.
    double px = double(probe.x), py = double(probe.y), pz = double(probe.z);

    WNV wnv = point_in_meshes_robust(px, py, pz, polygons, num_meshes);

    if (probe_side > 0)
        return wnv;
    else
    {
        for (size_t k = 0; k < wnv.size() && k < host_delta_w.size(); k++)
            wnv[k] -= host_delta_w[k];
        return wnv;
    }
}

// For reference propagation during subdivision
inline WNV classify_probe_by_raycast(
    pos_t const& probe, int /*ray_axis*/, bool /*ray_positive*/,
    IAABB const& /*aabb*/,
    std::vector<ConvexPolygon> const& polygons,
    int num_meshes)
{
    double px = double(probe.x), py = double(probe.y), pz = double(probe.z);
    return point_in_meshes_robust(px, py, pz, polygons, num_meshes);
}

} // namespace ember
